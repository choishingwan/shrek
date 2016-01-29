#include "genotypefilehandler.h"

GenotypeFileHandler::GenotypeFileHandler(){}
GenotypeFileHandler::~GenotypeFileHandler(){}


void GenotypeFileHandler::initialize(const Command &commander, const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Interval> &blockInfo){
    m_thread = commander.getThread();
    m_outPrefix = commander.getOutputPrefix();
    m_genotypeFilePrefix = commander.getLdFilePrefix();
    bool validate = commander.validate();
    bool mafFilt = commander.mafFilter();
    double mafThreshold = commander.getMaf();
    std::string line;
    if(commander.conRisk() || commander.diRisk()){
        std::string genotypeBim = commander.getGenotype()+".bim";
        std::string genotypeBed = commander.getGenotype()+".bed";
        bool bfile_SNP_major = openPlinkBinaryFile(genotypeBed, m_genoFile);
        if(bfile_SNP_major){
        //This is ok
        }
        else{
            throw std::runtime_error("We currently have no plan of implementing the individual-major mode. Please use the snp-major format");
        }

    }
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
    size_t warnings=0;
    size_t mafFilteredSnp=0;
    size_t finalNumSnp=0;
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
                    snpLoc = snpIndex.at(rs);
                    if(!validate || snpList.at(snpLoc).Concordant(chr, bp, rs)){
                        if(!validate && !snpList.at(snpLoc).Concordant(chr, bp, rs)) warnings++;
                        snpList.at(snpLoc).setFlag(0, true); //Now that it is in the genotype file, it contains the LD info.
                        m_inclusion.back()=snpLoc;
                        duplicateCheck[rs] = true;
                        finalNumSnp++;
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
                if(mafFilt &&  mafThreshold > currentMaf && snp){
                    m_inclusion.back()= -1;
                    mafFilteredSnp++;
                    finalNumSnp--;
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
    std::cerr << std::endl;
    std::cerr << "Linkage File information: " << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "Number of Samples: " << m_ldSampleSize << std::endl;
    std::cerr << "Duplicated SNPs:   " << duplicateCount << std::endl;
    std::cerr << "Invalid SNPs:      " << ignoreSnp << std::endl;
    std::cerr << "Filtered SNPs:     " << mafFilteredSnp << std::endl;
    std::cerr << "Final SNPs number: " << finalNumSnp << std::endl;
    if(warnings!=0) std::cerr << "WARNING: " << warnings << " invalid SNPs included" << std::endl;
    //Now we need to use the m_inclusion vector and the bim file to get get the intervals
    buildBlocks(bimFileName, blockInfo, commander.getDistance());
    m_finalSnpNumber =finalNumSnp;
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

    /**
     *  Concept of blockInfo:
     *  Start = the current block start (inclusive)
     *  End = the current block end (exclusive)
     *  So for the last block, the block end = size of SNP;
     */
    /**
     *  We only care about SNPs that were included.
     */
    std::string currentChr ="";
    std::string line;
    size_t currentIndex = 0;
    size_t prevIndex = 0;
    size_t startIndex = 0; //The beginning index of the current block
    size_t startBp = 0; //The beginning bp of the current block
    size_t prevBp = 0; //The bp location of the previous SNP
    bool started = false;


    while(std::getline(bimFile, line)){
        line = usefulTools::trim(line);
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, "\t ", &token);
            if(token.size() >= 6 && m_inclusion[currentIndex]!=-1){
                //Only work when this SNP is required
                std::string chr= token[0];
                size_t bp = std::atoi(token[3].c_str());
                if(!started){
                    //First, so just add everything
                    currentChr = chr;
                    startBp = bp;
                    startIndex = currentIndex;
                    started = true;
                }
                else if(chr.compare(currentChr)!=0){
                    //This is a new chromosome
                    //Check if there is any stuff here first
                    if(blockInfo.empty()){
                        //We need to use the previous information to build the first block
                        blockInfo.push_back(new Interval(currentChr, startIndex, prevIndex));
                    }
                    else if(blockInfo[blockInfo.size()-1].getChr().compare(currentChr)!= 0){
                        //That was a very small chromosome, still need to build it
                        blockInfo.push_back(new Interval(currentChr, startIndex, prevIndex));
                    }
                    else{
                        //We can directly extend the last block from the previous chromosome
                        blockInfo[blockInfo.size()-1].setEnd(prevIndex);
                    }
                    currentChr = chr;
                    startBp = bp; //Update the current location of the SNP
                    startIndex = currentIndex;
                }
                else if(bp-startBp > distance){
                    //They are on the same chromosome and there are two conditions
                    if(bp-prevBp > distance){
                        //Way too far away, basically they will immediately form a block together
                        //In this case, we want to extend the last block instead of building a new block
                        if(blockInfo.empty() || blockInfo[blockInfo.size()-1].getChr().compare(currentChr)!= 0){
                            //This is the first block of this chromosome
                            //We use prevIndex as the end here such that the end of this block != the start of next block, indicating they were too far
                            blockInfo.push_back(new Interval(currentChr, startIndex, prevIndex));
                        }
                        else{
                            //We will directly extend the previous interval
                            blockInfo[blockInfo.size()-1].setEnd(prevIndex);
                        }
                        startIndex = currentIndex;
                        startBp = bp;

                    }
                    else{
                        //This is the normal situation
                        //Because they are still continuous, we use currentIndex as the end where the end of this block = the start of next block
                        blockInfo.push_back(new Interval(currentChr, startIndex, currentIndex));
                        //So the block should start at the last SNP start and end with the previous SNP
                        startIndex = currentIndex;
                        startBp = bp;
                    }
                }
                currentIndex++;
                prevBp = bp;
                prevIndex = currentIndex; //Because it is exclusion, so it will be +1 already.
            }
            else if(token.size()>=6){
                currentIndex++; //This is the location of the SNP in the file, will increment no matter what happened
            }
        }

    }
    bimFile.close();
    if(prevIndex != blockInfo[blockInfo.size()-1].getEnd()){
        //It basically means that we have finished the whole file and still have SNPs not assigned to blocks
        //Only situation where this is possible is when bp - startBp < distance
        if(blockInfo.empty() || blockInfo[blockInfo.size()-1].getChr().compare(currentChr)!= 0){
            //the second condition is unlikely, but we will put it there for safety
            blockInfo.push_back(new Interval(currentChr, startIndex, prevIndex));
        }
        else {
            blockInfo[blockInfo.size()-1].setEnd(prevIndex);
        }
    }
//
//    for(size_t i = 0; i < blockInfo.size(); ++i){
//        std::cerr << blockInfo[i].getChr() << "\t" << blockInfo[i].getStart() << "\t" << blockInfo[i].getEnd() << std::endl;
//    }



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
	else if ( !v1_bfile ){
		if ( b[0] ) bfile_SNP_major = true;
		else bfile_SNP_major = false;
		std::cerr << "Binary PED file is v0.99" << std::endl;
		if (bfile_SNP_major) std::cerr << "Detected that binary PED file is in SNP-major mode" << std::endl;
		else std::cerr << "Detected that binary PED file is in individual-major mode" << std::endl;
	}
	return bfile_SNP_major;
}



void GenotypeFileHandler::getSnps(boost::ptr_vector<Genotype> &genotype, std::deque<size_t> &snpLoc, std::deque<size_t> &ldLoc, bool &chromosomeStart, bool &chromosomeEnd, size_t &prevResidual, boost::ptr_vector<Interval> &blockInfo, const std::vector<int> &genoInclusion, Eigen::MatrixXd &betaInfo){
    size_t startRange = blockInfo[prevResidual].getStart();
    size_t endRange=0;
    std::string currentChr = blockInfo[prevResidual].getChr();
    size_t i = prevResidual;
    size_t range = m_thread;
    if(chromosomeStart)range +=2;
    size_t loopEnd = prevResidual+range;
    bool firstEntry = true;
    for(; i < loopEnd && i < blockInfo.size(); ++i){
        if(blockInfo[i].getChr().compare(currentChr)!=0){
            //We are heading for another chromosome
            chromosomeEnd = true;
            break;
        }else if(!firstEntry && blockInfo[i].getStart() != endRange){
            std::cerr << "Not equal" << std::endl;
            firstEntry = false;
            //Again, we are working on another independent fragment
            chromosomeEnd = true; //Not really a new chromosome, but the behaviour of the programme should be the same
            break;
        }
        else{
            endRange = blockInfo[i].getEnd();
        }
    }
    if(i == blockInfo.size() || (i < blockInfo.size() && blockInfo[i].getChr().compare(currentChr)!=0)){ //The next block is from another chromosome
        chromosomeEnd=true;
    }
    prevResidual=i;
    /** This whole thing should basically be the same, except we now need to get additional genotypes **/
    /** We do assume the sort order for both the LD file and the geno file to be the same */
    //So we will get all the required SNPs within this region
    //for each required SNPs, we will use m_genoIter to iterate to the location.



    /** We use m_genoIter to indicate the index of the current genotype **/
    /** IMPORTANT: genotype is actually for storing the genotypes for LD calculation, not the sample genotype
     *             the sample genotype should be stored in betaInfo instead
     */
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
		m_snpIter++;
    }



}

void GenotypeFileHandler::getSnps(boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::deque<size_t> &ldLoc, bool &chromosomeStart, bool &chromosomeEnd, size_t &prevResidual, boost::ptr_vector<Interval> &blockInfo){
    //If this is the start of chromosome, we need to process one more bin
    //Otherwise we will process thread bins
    size_t startRange = blockInfo[prevResidual].getStart();
    size_t endRange=0;
    std::string currentChr = blockInfo[prevResidual].getChr();
    size_t i = prevResidual;
    size_t range = m_thread;
    if(chromosomeStart)range +=2;

    size_t loopEnd = prevResidual+range;
    /**
     *  Important property here is that in this for loop, i will +1 if it is a normal exit
     *  but will remain unchanged if it was break;
     */
    bool firstEntry = true;
    for(; i < loopEnd && i < blockInfo.size(); ++i){
        if(blockInfo[i].getChr().compare(currentChr)!=0){
            //We are heading for another chromosome
            chromosomeEnd = true;
            break;
        }else if(!firstEntry && blockInfo[i].getStart() != endRange){
            std::cerr << "Not equal" << std::endl;
            firstEntry = false;
            //Again, we are working on another independent fragment
            chromosomeEnd = true; //Not really a new chromosome, but the behaviour of the programme should be the same
            break;
        }
        else{
            endRange = blockInfo[i].getEnd();
        }
    }
    //i will be +1 because of the nature of for loop
    if(i == blockInfo.size() || (i < blockInfo.size() && blockInfo[i].getChr().compare(currentChr)!=0)){ //The next block is from another chromosome
        chromosomeEnd=true;
    }
    /**
     *  As a result of property of for loop, the prevResidual will ALWAYS be pointing to the next block.
     */
    prevResidual=i;
    //So we will get all the required SNPs within this region
    while(m_snpIter < m_inputSnp && m_snpIter <= endRange){
            //std::cerr << "reading: " << m_snpIter << std::endl;
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
		m_snpIter++;
    }
}
