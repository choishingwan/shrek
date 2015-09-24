#include "genotypefilehandler.h"

void GenotypeFileHandler::initialize(Command *commander, const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Interval> &blockInfo){
    m_thread = commander->getThread();
    m_outPrefix = commander->getOutputPrefix();
    m_genotypeFilePrefix = commander->getLdFilePrefix();
    bool validate = commander->validate();
    bool mafFilt = commander->mafFilter();
    double mafThreshold = commander->getMaf();
    std::string line;
    //Get the number of samples in the ld file
    std::string famFileName = m_genotypeFilePrefix+".fam";
    std::ifstream famFile;
    famFile.open(famFileName.c_str());
    if(!famFile.is_open()){
        throw std::runtime_error("Cannot open fam file");
    }
    while(std::getline(famFile, line)){
        line = usefulTools::trim(line);
        if(!line.empty()) m_ldSampleSize++;
    }
    famFile.close();
    //Check the bed file, and filter all SNPs not passing the MAF filtering
    std::string bedFileName = m_genotypeFilePrefix+".bed";
	bool bfile_SNP_major = openPlinkBinaryFile(bedFileName, m_bedFile); //We will try to open the connection to bedFile
    if(bfile_SNP_major){
        //This is ok
    }
    else{
        throw std::runtime_error("We currently have no plan of implementing the individual-major mode. Please use the snp-major format");
    }
    //We need to know the number of SNPs when transversing the bed file
    std::string bimFileName = m_genotypeFilePrefix+".bim";
    std::ifstream bimFile;
    bimFile.open(bimFileName.c_str());
    if(!bimFile.is_open()){
        throw std::runtime_error("Cannot open bim file");
    }
	std::string prevChr="";
	std::map<std::string, bool> duplicateCheck,sortCheck;
    size_t duplicateCount= 0;
    size_t ignoreSnp=0;
    size_t mafFilteredSnp=0;

    while(std::getline(bimFile,line)){
        line = usefulTools::trim(line);
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, "\t ", &token);
            if(token.size() >=6){
                std::string chr= token[0];
                std::string rs = token[1];
                //if(m_chrCount.find(chr)==m_chrCount.end()) m_chrCount[chr]=1;
                //else m_chrCount[chr]++;
                size_t bp = std::atoi(token[3].c_str());
                m_inputSnp++; //Number of total input Snps
                int snpLoc =-1;
                m_inclusion.push_back(-1); //Default is not including the snp
                //First run check a few things,
                //1. the bim is sorted
                //2. indicate which SNPs are required.
                if(prevChr.empty()){
                    sortCheck[chr]=true;
                }
                else if(prevChr.compare(chr) != 0){
                    if(sortCheck.find(chr)!=sortCheck.end()){
                        throw std::runtime_error("The programme require the SNPs to be sorted according to their chromosome.");
                    }
                    else sortCheck[chr] = true;
                }

                if(snpIndex.find(rs)!=snpIndex.end() && duplicateCheck.find(rs)==duplicateCheck.end()){
                    //This is something that we need
                    if(!validate || snpList.at(snpLoc).Concordant(chr, bp, rs)){
                        snpLoc = snpIndex.at(rs);
                        snpList.at(snpLoc).setFlag(0, true); //Now that it is in the genotype file, it contains the LD info.
                        m_inclusion.back()=snpLoc;
                        duplicateCheck[rs] = true;
                    }
                    else{
                        ignoreSnp++;
                    }
                }
                else if(duplicateCheck.find(rs) != duplicateCheck.end()){
                    duplicateCount++;
                }
                //Now perform the MAF check
                bool snp = false;
                if(m_inclusion.back() != -1){//indicate whether if we need this snp
                    snp=true;
                }
                size_t indx = 0; //The iterative count
                size_t alleleCount=0;
                while ( indx < m_ldSampleSize ){
                    std::bitset<8> b; //Initiate the bit array
                    char ch[1];
                    m_bedFile.read(ch,1); //Read the information
                    if (!m_bedFile){
                        throw "Problem with the BED file...has the FAM/BIM file been changed?";
                    }
                    b = ch[0];
                    int c=0;
                    while (c<7 && indx < m_ldSampleSize ){ //Going through the bit flag. Stop when it have read all the samples as the end == NULL
				//As each bit flag can only have 8 numbers, we need to move to the next bit flag to continue
                        ++indx; //so that we only need to modify the indx when adding samples but not in the mean and variance calculation
                        if (snp){
                            int first = b[c++];
                            int second = b[c++];
                            if(first == 1 && second == 0) first = 0; //We consider the missing value to be reference
                            alleleCount += first+second;
                        }
                        else{
                            c+=2;
                        }
                    }
                }
                double currentMaf = (alleleCount+0.0)/(2.0*m_ldSampleSize*1.0);
                currentMaf = (currentMaf > 0.5)? 1-currentMaf : currentMaf;
                if(mafFilt &&  mafThreshold > currentMaf){
                    m_inclusion.back()= -1;
                    mafFilteredSnp++;
                }
                //now the m_inclusion should include all the information we need: The SNP location if they are included and what snps to ignore
            }
            else{
                throw std::runtime_error("Malformed bim file, please check your input");
            }
        }
    }
    bimFile.close();
    m_bedFile.close();
    std::cerr << "Linkage File information: " << std::endl;
    std::cerr << "Duplicated SNPs: " << duplicateCount << std::endl;
    std::cerr << "Invalid SNPs:    " << ignoreSnp << std::endl;
    std::cerr << "Filtered SNPs:   " << mafFilteredSnp << std::endl;
    //Now we need to use the m_inclusion vector and the bim file to get get the intervals
    buildBlocks(bimFileName, blockInfo, commander->getDistance());
    //Now open the bed file to prepare for whatever happen next
	bfile_SNP_major = openPlinkBinaryFile(bedFileName, m_bedFile); //We will try to open the connection to bedFile
    if(bfile_SNP_major){
        //This is ok
    }
    else{
        throw std::runtime_error("We currently have no plan of implementing the individual-major mode. Please use the snp-major format");
    }

}

void GenotypeFileHandler::buildBlocks(std::string bimFileName, boost::ptr_vector<Interval> &blockInfo, size_t distance){
    //Should build the block here
    //We will merge the last block such that the last block will still be > distance
    std::ifstream bimFile;
    bimFile.open(bimFileName.c_str());
    if(!bimFile.is_open()){
        throw std::runtime_error("Cannot open bim file");
    }
    std::string currentChr ="";
    std::string line;
    size_t currentIndex = 0;
    size_t currentStart = 0;
    size_t prevLoc = 0;
    bool started = false;
    while(std::getline(bimFile, line)){
        line = usefulTools::trim(line);
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, "\t ", &token);
            //The m_inclusion will be -1 if we don't need the SNP
            if(token.size() >= 6 && m_inclusion[currentIndex]!=-1){
                std::string chr= token[0];
                size_t bp = std::atoi(token[3].c_str());
                if(!started){
                    currentChr = chr;
                    currentStart = bp;
                    started = true;
                }
                else if(chr.compare(currentChr)!= 0){
                    //new chromosome
                    blockInfo[blockInfo.size()-1].setEnd(prevLoc);
                    //blockInfo.push_back(new Interval(currentChr, currentStart, prevLoc));
                    currentChr = chr;
                    currentStart = bp;
                }
                else if(bp-currentStart > distance){
                    //Now we are off, so this should be the interval
                    blockInfo.push_back(new Interval(currentChr, currentStart, prevLoc));
                    currentStart = prevLoc;
                    if(bp-currentStart > distance){
                        //The new SNP is also <distance> away from the last SNP
                        blockInfo.push_back(new Interval(currentChr, currentStart, bp));
                        currentStart = bp;
                    }
                }
                prevLoc = bp;
                currentIndex++;
            }
        }
    }
    bimFile.close();
    if(prevLoc != blockInfo[blockInfo.size()-1].getEnd()){
        blockInfo[blockInfo.size()-1].setEnd(prevLoc);
    }
}

bool GenotypeFileHandler::openPlinkBinaryFile(const std::string s, std::ifstream & BIT){
	BIT.open(s.c_str(), std::ios::in | std::ios::binary);
	if(!BIT.is_open()){
        throw "Cannot open the bed file";
	}
	//std::cerr << "BIT open" << std::endl;
	// 1) Check for magic number
	// 2) else check for 0.99 SNP/Ind coding
	// 3) else print warning that file is too old
	char ch[1];
	BIT.read(ch,1);
	std::bitset<8> b;
	b = ch[0];
	bool bfile_SNP_major = false;
	bool v1_bfile = true;
	// If v1.00 file format
	// Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
	//std::cerr << "check magic number" << std::endl;
	if (   ( b[2] && b[3] && b[5] && b[6] ) &&
       ! ( b[0] || b[1] || b[4] || b[7] )    ){
	// Next number
	BIT.read(ch,1);
	b = ch[0];
	if (   ( b[0] && b[1] && b[3] && b[4] ) &&
          ! ( b[2] || b[5] || b[6] || b[7] )    ){
			// Read SNP/Ind major coding
			BIT.read(ch,1);
			b = ch[0];
			if ( b[0] ) bfile_SNP_major = true;
			else bfile_SNP_major = false;

			//if (bfile_SNP_major) std::cerr << "Detected that binary PED file is v1.00 SNP-major mode" << std::endl;
			//else std::cerr << "Detected that binary PED file is v1.00 individual-major mode" << std::endl;

		} else v1_bfile = false;

	} else v1_bfile = false;
	// Reset file if < v1
	if ( ! v1_bfile ) {
		std::cerr << "Warning, old BED file <v1.00 : will try to recover..." << std::endl;
		std::cerr << "  but you should --make-bed from PED )" << std::endl;
		BIT.close();
		BIT.clear();
		BIT.open(s.c_str(), std::ios::in | std::ios::binary);
		BIT.read(ch,1);
		b = ch[0];
	}
	// If 0.99 file format
	if ( (!v1_bfile) && ( b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7] ) ){
		std::cerr << std::endl << " *** Possible problem: guessing that BED is < v0.99      *** " << std::endl;
		std::cerr << " *** High chance of data corruption, spurious results    *** " << std::endl;
		std::cerr << " *** Unless you are _sure_ this really is an old BED file *** " << std::endl;
		std::cerr << " *** you should recreate PED -> BED                      *** " << std::endl << std::endl;
		bfile_SNP_major = false;
		BIT.close();
		BIT.clear();
		BIT.open(s.c_str(), std::ios::in | std::ios::binary);
	}
	else if ( ! v1_bfile ){
		if ( b[0] ) bfile_SNP_major = true;
		else bfile_SNP_major = false;
		std::cerr << "Binary PED file is v0.99" << std::endl;
		if (bfile_SNP_major) std::cerr << "Detected that binary PED file is in SNP-major mode" << std::endl;
		else std::cerr << "Detected that binary PED file is in individual-major mode" << std::endl;
	}
	return bfile_SNP_major;
}

void GenotypeFileHandler::getSnps(boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::deque<size_t> &ldLoc, boost::ptr_vector<Snp> &snpList, bool &chromosomeStart, bool &chromosomeEnd, size_t &prevResidual, boost::ptr_vector<Interval> &blockInfo){
    //If this is the start of chromosome, we need to process one more bin
    //Otherwise we will process thread bins
    size_t startRange = blockInfo[prevResidual].getStart();
    size_t endRange=0;
    std::string currentChr = blockInfo[prevResidual].getChr();
    size_t i = prevResidual;
    size_t range = m_thread;
    if(chromosomeStart)range +=2;
    for(; i< prevResidual+range && i < blockInfo.size(); ++i){
            if(blockInfo[i].getChr().compare(currentChr)!=0){
            i= i-1; //This mean we working on the last block of this chromosome
            chromosomeEnd = true;
            break;
        }
    }
    endRange = blockInfo[i].getEnd();
    //So we will get all the required SNPs within this region

    while(m_snpIter < m_inputSnp && m_snpIter <= endRange){
        bool snp = false;
		if(m_inclusion[m_snpIter] != -1){//indicate whether if we need this snp
            if(m_snpIter < startRange){
                std::logic_error("Warning, something went wrong with my logic, really sorry");
            }
			snp=true;
		}
        if(snp){
            genotype.push_back(new Genotype());
            snpLoc.push_back(m_inclusion[m_snpIter]); //need better way, such that we know the block problem
            ldLoc.push_back(m_snpIter);
		}
        size_t indx = 0; //The iterative count
        double oldM=0.0, newM=0.0,oldS=0.0, newS=0.0;
        size_t alleleCount=0;
		while ( indx < m_ldSampleSize ){
			std::bitset<8> b; //Initiate the bit array
			char ch[1];
			m_bedFile.read(ch,1); //Read the information
			if (!m_bedFile){
				throw std::runtime_error("Problem with the BED file...has the FAM/BIM file been changed?");
			}
			b = ch[0];
			int c=0;
			while (c<7 && indx < m_ldSampleSize ){ //Going through the bit flag. Stop when it have read all the samples as the end == NULL
				//As each bit flag can only have 8 numbers, we need to move to the next bit flag to continue
                ++indx; //so that we only need to modify the indx when adding samples but not in the mean and variance calculation
				if (snp){
					int first = b[c++];
					int second = b[c++];
					if(first == 1 && second == 0) first = 3; //Missing value should be 3
					genotype.back().AddsampleGenotype(first+second, indx-1); //0 1 2 or 3 where 3 is missing
					alleleCount += first+second;
					double value = first+second+0.0;
                    if(indx==1){
                        oldM = newM = value;
                        oldS = 0.0;
                    }
                    else{
                        newM = oldM + (value-oldM)/(indx);
                        newS = oldS + (value-oldM)*(value-newM);
                        oldM = newM;
                        oldS = newS;
                    }

				}
				else{
					c+=2;
				}
			}
		}
		if(snp){
			indx > 0 ? genotype.back().Setmean(newM) : genotype.back().Setmean(0.0);
			indx > 1 ? genotype.back().SetstandardDeviation(std::sqrt(newS/(indx - 1.0))) : genotype.back().SetstandardDeviation(0.0);
		}
    }
}
