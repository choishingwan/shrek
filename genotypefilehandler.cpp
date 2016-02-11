#include "genotypefilehandler.h"

GenotypeFileHandler::GenotypeFileHandler()
{
    //ctor
}

GenotypeFileHandler::~GenotypeFileHandler()
{
    //dtor
    if(m_buffGenotype!=nullptr)
        delete m_buffGenotype;
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





void GenotypeFileHandler::initialize(const Command &commander, const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList){
    m_keepAmbiguous = commander.keepAmbiguous();
    m_thread = commander.getNThread();
    m_genotypeFilePrefix = commander.getReferenceFilePrefix();
    m_mafThreshold = commander.getMAF();
    m_blockSize = commander.getSizeOfBlock();
    // First, we need to know the number of samples such that we can read the plink binary file correctly
    std::string famFileName = m_genotypeFilePrefix+".fam";
    std::ifstream famFile;
    famFile.open(famFileName.c_str());
    if(!famFile.is_open()){
        throw std::runtime_error("Cannot open fam file");
    }
    std::string line;
    while(std::getline(famFile, line)){
        line = usefulTools::trim(line);
        if(!line.empty()) m_nRefSample++;
    }
    famFile.close();
    // We need to tell Genotype what is the maximum number of samples that can present in the reference panel
    // therefore we must do
    Genotype::SetsampleNum(m_nRefSample);
    // End up realized that if I want to have a dynamic block construction
    // when I read the SNPs, I will still need to read both the bim and bed
    // files together, so don't bother with reading the file now

    // However, I should still initialize the bed file

    std::string bedFileName = m_genotypeFilePrefix+".bed";

    bool bfile_SNP_major = openPlinkBinaryFile(bedFileName, m_bedFile); //We will try to open the connection to bedFile
    if(!bfile_SNP_major)
        throw std::runtime_error("We currently have no plan of implementing the individual-major mode. Please use the snp-major format");
    std::string bimFileName = m_genotypeFilePrefix+".bim";
    std::ifstream bimFile;
    bimFile.open(bimFileName.c_str());
    if(!bimFile.is_open())  throw std::runtime_error("Cannot open bim file");
//    m_bimFile.open(bimFileName.c_str());
//    if(!m_bimFile.is_open()){
//        throw std::runtime_error("Cannot open bim file");
//    }
//    fprintf(stderr, "Start reading the first SNP\n");
//    initializeSNP(snpIndex, snpList);
//    fprintf(stderr, "Finished initializing the first SNP\n");

    size_t emptyLineCheck = 0; // Only use for checking for empty line, consider this as a safeguard
    while(std::getline(bimFile, line)){
        line = usefulTools::trim(line);
        //The only thing I cannot filter beforehand is the MAF
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, "\t", &token);
            m_inclusion.push_back(-1);
            if(token.size()>=6){
                std::string rs = token[1];
                if(m_duplicateCheck.find(rs)!=m_duplicateCheck.end())   m_nDuplicated++; // Don't want it if it is duplicated (impossible here actually)
                else if(snpIndex.find(rs)!=snpIndex.end()){
                    size_t location = snpIndex.at(rs);
                    bool ambig = false;
                    if(snpList.at(location).concordant(token[0], atoi(token[3].c_str()), rs, token[4], token[5], ambig )){
                        if(!m_keepAmbiguous && ambig) m_nAmbig ++; // We don't want to keep ambiguous SNPs and this is an ambiguous SNP
                        else m_inclusion.back()=location;
                    }
                    else m_nInvalid++;
                }
            }
            else throw "Malformed bim file, does not contain all necessary information";
        }
        else emptyLineCheck++;
    }
    bimFile.close();
    fprintf(stderr,"\n Reference Panel Information:\n");
    fprintf(stderr,"=================================\n");
    fprintf(stderr, "%lu samples were found in the reference panel\n", m_nRefSample);
    fprintf(stderr, "%lu discordant SNP(s) identified\n", m_nInvalid);
    if(m_nAmbig!=0)
        fprintf(stderr, "%lu ambiguous SNP(s) filtered\n", m_nAmbig);
    if(m_nDuplicated)
        fprintf(stderr, "%lu duplicated SNP(s) found in reference panel\n", m_nInvalid);
    m_snpIter=0;
    m_nSnp = m_inclusion.size(); //reduce the number of time required for accessing the size
}



// Find the first SNP that for us to include.
// Then read all the SNPs within the region

void GenotypeFileHandler::getBlock(boost::ptr_vector<Snp> &snpList, boost::ptr_list<Genotype> &genotype, std::list<size_t> &snpLoc, bool &finalizeBuff, bool &completed, std::deque<std::list<size_t>::iterator > &boundary){
    // We try to use the old version of stuff to work on
    // Start of block info


    bool starting  = true;
    int blockStartIndex = 0;
    std::string blockChr="";
    size_t blockStartLoc=0;
    size_t lastUsedLoc =0;
    for(; m_snpIter < m_nSnp; ++m_snpIter){
        bool snp = false;
        if(starting){
            blockStartIndex = m_inclusion[m_snpIter];
            if(blockStartIndex!= -1){
                blockChr = snpList.at(blockStartIndex).getChr();
                blockStartLoc = snpList.at(blockStartIndex).getLoc();
                lastUsedLoc = blockStartLoc;
                snp = true;
            }
        }
        //indicate whether if we need this snp
		if(m_inclusion[m_snpIter] != -1) snp=true;
		if(snp){
            size_t currentLoc = snpList.at(m_inclusion[m_snpIter]).getLoc();

            std::string currentChr = snpList.at(m_inclusion[m_snpIter]).getChr();
            if(currentChr.compare(blockChr)!=0 ||  // new chromosome
                currentLoc - blockStartLoc > m_blockSize || // way too far
                currentLoc - lastUsedLoc > m_blockSize){
                // we don't need this SNP, so we will return immediately
                finalizeBuff=true;
                return;
            }
            else{
                // This is something we need
                Genotype *tempGenotype = new Genotype();
                size_t indx = 0; //The iterative count
                double oldM=0.0, newM=0.0,oldS=0.0, newS=0.0;
                size_t alleleCount=0, validSample=0;
                while ( indx < m_nRefSample ){
                    std::bitset<8> b; //Initiate the bit array
                    char ch[1];
                    m_bedFile.read(ch,1); //Read the information
                    if (!m_bedFile) throw std::runtime_error("Problem with the BED file...has the FAM/BIM file been changed?");
                    b = ch[0];
                    int c=0;
                    while (c<7 && indx < m_nRefSample ){ //Going through the bit flag. Stop when it have read all the samples as the end == NULL
                        //As each bit flag can only have 8 numbers, we need to move to the next bit flag to continue
                        ++indx; //so that we only need to modify the indx when adding samples but not in the mean and variance calculation
                        int first = b[c++];
                        int second = b[c++];

                        if(first == 1 && second == 0){
                            first = 3; //Missing value should be 3
                        }
                        else{
                            validSample++;
                            alleleCount += first+second;
                            double value = first+second+0.0;
                            if(validSample==1){
                                oldM = newM = value;
                                oldS = 0.0;
                            }
                            else{
                                newM = oldM + (value-oldM)/(validSample);
                                newS = oldS + (value-oldM)*(value-newM);
                                oldM = newM;
                                oldS = newS;
                            }
                        }
                        tempGenotype->AddsampleGenotype(first+second, indx-1); //0 1 2 or 3 where 3 is missing
                    }
                }
                validSample > 0 ? tempGenotype->Setmean(newM) : tempGenotype->Setmean(0.0);
                validSample > 1 ? tempGenotype->SetstandardDeviation(std::sqrt(newS/(validSample - 1.0))) : tempGenotype->SetstandardDeviation(0.0);
                double currentMaf = (alleleCount+0.0)/(2.0*validSample*1.0);
                currentMaf = (currentMaf > 0.5)? 1-currentMaf : currentMaf;
                if(m_mafThreshold > currentMaf){
                    m_inclusion[m_snpIter]=-1;
                    m_nFilter++;
                    delete tempGenotype;
                }
                else{
                    genotype.push_back(tempGenotype);
                    snpLoc.push_back(currentLoc);
                    lastUsedLoc  = currentLoc;
                    // fprintf(stderr, "Check %lu\n",m_snpLoc);
                    if(starting) boundary.push_back(std::prev(snpLoc.end()));
                    starting = false;
                }
            }
		}
		else{
            //Doesn't matter, just read the thing and proceed
            size_t indx = 0; //The iterative count
            while ( indx < m_nRefSample ){
                std::bitset<8> b; //Initiate the bit array
                char ch[1];
                m_bedFile.read(ch,1); //Read the information
                if (!m_bedFile) throw std::runtime_error("Problem with the BED file...has the FAM/BIM file been changed?");
                b = ch[0];
                int c=0;
                while (c<7 && indx < m_nRefSample ){ //Going through the bit flag. Stop when it have read all the samples as the end == NULL
                        //As each bit flag can only have 8 numbers, we need to move to the next bit flag to continue
                    ++indx; //so that we only need to modify the indx when adding samples but not in the mean and variance calculation
                    c+=2;
                }
            }
		}
    }

    finalizeBuff = true;
    completed = true;

}

void GenotypeFileHandler::getSNP(boost::ptr_vector<Snp> &snpList, boost::ptr_list<Genotype> &genotype, std::list<size_t> &snpLoc, bool &finalizeBuff, bool &completed, std::deque<std::list<size_t>::iterator > &boundary){
    // boundary size should at most be 4

    // Very complicated need to write carefully
    // Something we know
    // When we reached the end of the chromosome or when the distance between the two SNPs are too large
    // we will set finalizeBuff to true such that other function will start cleaning up
//    /** The content of the boundary is important so it is worth the time to document it
//     *  The boundary should contain the boundaries of the block. The start of the first
//     *  block will always be 0. Let the start of the second block be x, then the first
//     *  block will be [0,x) and second block will be [x,y). Also, the distance of
//     *  SNP[x] - SNP[0] >= defined distance
//     */
//    // If the boundary size is 0, then we need a lot of blocks (3), otherwise, we only
//    // need to read one additional block
//    //size_t nBlock = (boundary.size()==0)? 4:1; //we read one more, just so that we might need to merge the remaining information
//    /** From this point onward, we assume prevChr, prevLoc and bufferGenotype contains the last SNP entry **/
//    while( boundary.size() < 4){
//    //for(size_t i = 0; i < nBlock; ++i){
//        // Get get block here
//        // When trying to get the blocks, we will try to get all the SNPs within region of the SNP in the bufferGenotype,
//        // we will also read one extra SNP (e.g. first SNP off the bufferGenotype) and put it into the buffer.
//        getBlock(snpIndex, snpList, genotype, snpLoc, finalizeBuff,completed, boundary);
//        if(finalizeBuff && boundary.size() > 3){
//            // Now check if we need to merge the last two blocks
//            size_t lastSnpOfThirdBlock = (*boundary.back()); // The last of boundary indicate the start of the last block
//            lastSnpOfThirdBlock = snpList.at(lastSnpOfThirdBlock).getLoc();
//            size_t lastSnp = snpLoc.back(); //The last of snpLoc is the last of the last block
//            lastSnp = snpList.at(lastSnp).getLoc();
//            if(lastSnp-lastSnpOfThirdBlock <= m_blockSize){
//                boundary.pop_back();
//            }
//            return; //Done with this
//        }
//        else if(finalizeBuff){ // When we basically have too little SNPs as an input
//            while(boundary.size()!=1) boundary.pop_back();
//            return;
//        }
//    }
 /** I have updated the functions, now getSnp should just get more SNPs to fill in the last block **/


    size_t lastStartIndex = *(boundary.end());
    std::string blockChr = snpList.at(lastStartIndex).getChr();
    size_t blockStartLoc = snpList.at(lastStartIndex).getLoc();
    size_t lastUsedLoc = blockStartLoc;

    for(; m_snpIter < m_nSnp; ++m_snpIter){
        bool snp = false;
        //indicate whether if we need this snp
		if(m_inclusion[m_snpIter] != -1) snp=true;
		if(snp){
            size_t currentLoc = snpList.at(m_inclusion[m_snpIter]).getLoc();

            std::string currentChr = snpList.at(m_inclusion[m_snpIter]).getChr();
            if(currentChr.compare(blockChr)!=0 ||  // new chromosome
                currentLoc - blockStartLoc > m_blockSize || // way too far
                currentLoc - lastUsedLoc > m_blockSize){
                // we don't need this SNP, so we will return immediately
                finalizeBuff=true;
                return;
            }
            else{
                // This is something we need
                Genotype *tempGenotype = new Genotype();
                size_t indx = 0; //The iterative count
                double oldM=0.0, newM=0.0,oldS=0.0, newS=0.0;
                size_t alleleCount=0, validSample=0;
                while ( indx < m_nRefSample ){
                    std::bitset<8> b; //Initiate the bit array
                    char ch[1];
                    m_bedFile.read(ch,1); //Read the information
                    if (!m_bedFile) throw std::runtime_error("Problem with the BED file...has the FAM/BIM file been changed?");
                    b = ch[0];
                    int c=0;
                    while (c<7 && indx < m_nRefSample ){ //Going through the bit flag. Stop when it have read all the samples as the end == NULL
                        //As each bit flag can only have 8 numbers, we need to move to the next bit flag to continue
                        ++indx; //so that we only need to modify the indx when adding samples but not in the mean and variance calculation
                        int first = b[c++];
                        int second = b[c++];

                        if(first == 1 && second == 0){
                            first = 3; //Missing value should be 3
                        }
                        else{
                            validSample++;
                            alleleCount += first+second;
                            double value = first+second+0.0;
                            if(validSample==1){
                                oldM = newM = value;
                                oldS = 0.0;
                            }
                            else{
                                newM = oldM + (value-oldM)/(validSample);
                                newS = oldS + (value-oldM)*(value-newM);
                                oldM = newM;
                                oldS = newS;
                            }
                        }
                        tempGenotype->AddsampleGenotype(first+second, indx-1); //0 1 2 or 3 where 3 is missing
                    }
                }
                validSample > 0 ? tempGenotype->Setmean(newM) : tempGenotype->Setmean(0.0);
                validSample > 1 ? tempGenotype->SetstandardDeviation(std::sqrt(newS/(validSample - 1.0))) : tempGenotype->SetstandardDeviation(0.0);
                double currentMaf = (alleleCount+0.0)/(2.0*validSample*1.0);
                currentMaf = (currentMaf > 0.5)? 1-currentMaf : currentMaf;
                if(m_mafThreshold > currentMaf){
                    m_inclusion[m_snpIter]=-1;
                    m_nFilter++;
                    delete tempGenotype;
                }
                else{
                    genotype.push_back(tempGenotype);
                    snpLoc.push_back(currentLoc);
                    lastUsedLoc  = currentLoc;
                }
            }
		}
		else{
            //Doesn't matter, just read the thing and proceed
            size_t indx = 0; //The iterative count
            while ( indx < m_nRefSample ){
                std::bitset<8> b; //Initiate the bit array
                char ch[1];
                m_bedFile.read(ch,1); //Read the information
                if (!m_bedFile) throw std::runtime_error("Problem with the BED file...has the FAM/BIM file been changed?");
                b = ch[0];
                int c=0;
                while (c<7 && indx < m_nRefSample ){ //Going through the bit flag. Stop when it have read all the samples as the end == NULL
                        //As each bit flag can only have 8 numbers, we need to move to the next bit flag to continue
                    ++indx; //so that we only need to modify the indx when adding samples but not in the mean and variance calculation
                    c+=2;
                }
            }
		}
    }

    finalizeBuff = true;
    completed = true;
}


