#include "genotypefilehandler.h"

size_t GenotypeFileHandler::GetsampleSize() const { return m_ldSampleSize; }
size_t GenotypeFileHandler::GetestimateSnpTotal() const { return m_expectedNumberOfSnp; }

GenotypeFileHandler::GenotypeFileHandler(std::string genotypeFilePrefix, size_t thread, std::string outPrefix):m_genotypeFilePrefix(genotypeFilePrefix), m_thread(thread), m_outPrefix(outPrefix){
	m_defaultDistance=2000000;
	m_ldSampleSize = 0;
	m_expectedNumberOfSnp = 0;
    m_inputSnp =0;
	m_snpIter =0;
}


void GenotypeFileHandler::initialize(std::map<std::string, size_t> &snpIndex, std::vector<Snp*> *snpList, bool validate, bool maxBlockSet, size_t maxBlock, size_t minBlock, double const maf){
    size_t safeBlockRange = 1.0; //Use to multiply the #Snp in 1mb region to make sure the block will always include everything within the region
    std::string famFileName = m_genotypeFilePrefix +".fam";
    std::ifstream famFile;
    std::ofstream blockRecommendOut;
    bool stdOut = true;
    if(!m_outPrefix.empty()){
        stdOut = false;
        std::string blockRecName = m_outPrefix+".block";
        blockRecommendOut.open(blockRecName.c_str());
        if(!blockRecommendOut.is_open()){
            throw "Cannot open the block recommendation log file";
        }
    }
    famFile.open(famFileName.c_str());
    if(!famFile.is_open()){
        throw "Cannot open the fam file";
    }
    std::string line;
    while(std::getline(famFile, line)){
        line = usefulTools::trim(line);
        if(!line.empty()) m_ldSampleSize++;
    }
    famFile.close();
    std::cerr << "A total of " << m_ldSampleSize << " samples were found in the genotype file for LD construction" << std::endl << std::endl;
    Genotype::SetsampleNum(m_ldSampleSize);

    std::string bedFileName = m_genotypeFilePrefix+".bed";
	bool bfile_SNP_major = openPlinkBinaryFile(bedFileName, m_bedFile); //We will try to open the connection to bedFile
    if(bfile_SNP_major){
        //This is ok
    }
    else{
        throw "We currently have no plan of implementing the individual-major mode. Please use the snp-major format";
    }

    std::string bimFileName = m_genotypeFilePrefix+".bim";
    std::ifstream bimFile;
    bimFile.open(bimFileName.c_str());
    if(!bimFile.is_open()){
        throw "Cannot open bim file";
    }
	std::map<std::string, bool> duplicateCheck, sortCheck;
    int duplicateCount = 0;
    //check the optimum blockSize here

    if(stdOut) std::cerr << "chr\tRecommend\tFinal"<< std::endl;
    else blockRecommendOut << "chr\tRecommend\tFinal"<< std::endl;
    size_t currentMaxBlock = 0;
	std::deque<size_t> locList;
	std::string prevChr="";
	size_t prevSnpLoc = 0;
    while(std::getline(bimFile, line)){
        line = usefulTools::trim(line);
        if(!line.empty()){
            //m_expectedNumberOfSnp++;
            std::vector<std::string> token;
            usefulTools::tokenizer(line, "\t ", &token);
            if(token.size() >=6){
                std::string chr= token[0];
                std::string rs = token[1];
                if(m_chrCount.find(chr)==m_chrCount.end())m_chrCount[chr]=1;
                else m_chrCount[chr]++;
                size_t bp = std::atoi(token[3].c_str());
                m_inputSnp++; //Number of total input Snps
                int snpLoc =-1;
                m_inclusion.push_back(-1); //Default is not including the snp
                if(prevChr.empty()){
                    prevChr = chr;
                    prevSnpLoc = bp;
                    sortCheck[chr] = true;

                    m_chrExists.push_back(chr);
                    currentMaxBlock= 0;
                    if(snpIndex.find(rs)!=snpIndex.end()){
                        //This is something that we need
                        if(!validate || (*snpList).at(snpLoc)->Concordant(chr, bp, rs)){
                            snpLoc = snpIndex.at(rs);
                            (*snpList).at(snpLoc)->setFlag(0, true); //Now that it is in the genotype file, it contains the LD info.
                            m_inclusion.back()=snpLoc;
                            locList.push_back(bp);
                            duplicateCheck[rs] = true;
                        }
                        else{
                            std::cerr << rs << " has different information from p-value file! Will ignore Snp." << std::endl;
                            //basically we ignore the snp
                        }
                    }
                }
                else if(prevChr.compare(chr)==0){
                        if(bp < prevSnpLoc){
                            //The snps were not sorted
                            throw "The programme require the SNPs to be sorted according to their chromosome. Sorry.";
                        }
                        prevSnpLoc = bp;
                        if(snpIndex.find(rs)!=snpIndex.end() && duplicateCheck.find(rs)==duplicateCheck.end()){
                            if(!validate || (*snpList).at(snpLoc)->Concordant(chr, bp, rs)){
                                snpLoc = snpIndex.at(rs);
                                (*snpList).at(snpLoc)->setFlag(0, true); //Now that it is in the genotype file, it contains the LD info.
                                m_inclusion.back()=snpLoc;
                                if(bp-locList.front() > m_defaultDistance){
                                    if(currentMaxBlock < locList.size()*safeBlockRange) currentMaxBlock = locList.size()*safeBlockRange;
                                    while(bp-locList.front() > m_defaultDistance && !locList.empty()){
                                        locList.pop_front();
                                    }
                                }
                                locList.push_back(bp);
                            }
                            else{
                                std::cerr << rs << " has different information from p-value file! Will ignore Snp." << std::endl;
                                //basically we ignore the snp
                            }
                        }
                        else if(duplicateCheck.find(rs)!=duplicateCheck.end()){
                            duplicateCount++;
                        }

                }
                else if(prevChr.compare(chr) != 0){
                    if(sortCheck.find(chr)!=sortCheck.end()){
                        //Unsorted chromosome
                        throw "The programme require the SNPs to be sorted according to their chromosome. Sorry.";
                    }
                    //Valid new chromosome


                    if(currentMaxBlock < locList.size()*safeBlockRange) currentMaxBlock = locList.size()*safeBlockRange;
					if(currentMaxBlock%3 != 0){
                        currentMaxBlock = currentMaxBlock+3-currentMaxBlock%3;
					}
					if(stdOut) std::cerr <<prevChr << "\t" << currentMaxBlock <<"\t";
					else blockRecommendOut << prevChr << "\t" << currentMaxBlock <<"\t";
					if(maxBlockSet && maxBlock < currentMaxBlock ){
						currentMaxBlock = maxBlock;
					}
					if(currentMaxBlock < minBlock) currentMaxBlock = minBlock;
					if(stdOut) std::cerr << currentMaxBlock << std::endl;
					else blockRecommendOut << currentMaxBlock << std::endl;
					m_blockSizeTract[prevChr]=currentMaxBlock;
					locList.clear();
					//All new now
                    prevChr = chr;
                    prevSnpLoc = bp;
                    sortCheck[chr] = true;
                    m_chrExists.push_back(chr);
                    currentMaxBlock= 0;
                    if(snpIndex.find(rs)!=snpIndex.end()){
                        //This is something that we need
                        if(!validate || (*snpList).at(snpLoc)->Concordant(chr, bp, rs)){
                            snpLoc = snpIndex.at(rs);
                            (*snpList).at(snpLoc)->setFlag(0, true); //Now that it is in the genotype file, it contains the LD info.
                            m_inclusion.back()=snpLoc;
                            locList.push_back(bp);
                            duplicateCheck[rs] = true;
                        }
                        else{
                            std::cerr << rs << " has different information from p-value file! Will ignore Snp." << std::endl;
                            //basically we ignore the snp
                        }
                    }
                }
                //Try to do the maf check here

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

                double currentMaf = (alleleCount+0.0)/(2*m_ldSampleSize*1.0);
                currentMaf = (currentMaf > 0.5)? 1-currentMaf : currentMaf;
                //remove snps with maf too low
                if(maf >= 0.0 && maf > currentMaf){
                    m_inclusion.back()= -1;
                    std::cerr << "Snp: " << rs << " not included due to maf filtering" << std::endl;
                }
                else if(m_inclusion.back() != -1){
                    if(m_chrProcessCount.find(chr)==m_chrProcessCount.end()) m_chrProcessCount[chr] = 1;
                    else m_chrProcessCount[chr]++;
                }
            }
            else{
                throw "Line in bim file has incorrect number of dimension";
            }
        }
    }
    //Only do this if we still have Snps left
    if(currentMaxBlock < locList.size()*safeBlockRange) currentMaxBlock = locList.size()*safeBlockRange;
    if(currentMaxBlock%3 != 0){
        currentMaxBlock = currentMaxBlock+3-currentMaxBlock%3;
    }
    if(stdOut) std::cerr << prevChr << "\t" << currentMaxBlock << "\t";
    else blockRecommendOut << prevChr << "\t" << currentMaxBlock << "\t";
    if(maxBlockSet && maxBlock < currentMaxBlock ){
        currentMaxBlock = maxBlock;
    }
    if(currentMaxBlock < minBlock) currentMaxBlock = minBlock;
    if(stdOut) std::cerr << currentMaxBlock << std::endl;
    else blockRecommendOut << currentMaxBlock << std::endl;
    m_blockSizeTract[prevChr]=currentMaxBlock;
	//The reason of not using this check is because the currentMaxBlock is used for all chromosome
	//If the chromosome does not have anything in the p-value file, currentMaxBlock will most likely be 0
	//if(currentMaxBlock == 0) throw "The block size is 0, most likely your input file is problematic. Please check.";
    if(!stdOut) blockRecommendOut.close();

    bimFile.close();
    m_bedFile.close();
    if(duplicateCount == 0) std::cerr << "There are no duplicated snps in the LD file" << std::endl;
    else{
		std::cerr << "A total of " << duplicateCount << " Snps in the LD file were duplicated" << std::endl;
		std::cerr << "Only the first instance of each Snp will be used" << std::endl;
    }
	std::cerr << std::endl;
    bfile_SNP_major = openPlinkBinaryFile(bedFileName, m_bedFile); //We will try to open the connection to bedFile
    if(bfile_SNP_major){
        //This is ok
    }
    else{
        throw "We currently have no plan of implementing the individual-major mode. Please use the snp-major format";
    }
	//Initialize the variable for getSnps
	m_processed = 0;
	m_targetProcessed=0;
	//Now calculate the block start size and update the snps accordingly
    //For each chromosome we will need to have the block information

    //Check m_chrprocessCount
    size_t index = 0;
    for(size_t i = 0; i < m_chrExists.size(); ++i){
        std::string chr = m_chrExists[i];
        if(m_chrProcessCount.find(chr)!= m_chrProcessCount.end()){
            m_expectedNumberOfSnp +=m_chrProcessCount[chr];
            if(m_chrCount.find(chr)!=m_chrCount.end()){
                //Build the vector first
                std::vector<size_t> blockInfo(m_chrProcessCount[chr],3.0);
                size_t block  = m_blockSizeTract[chr]/3;
                if(block > blockInfo.size()){
                    std::fill_n(blockInfo.begin(), blockInfo.size(), 1.0);
                }
                else{
                    std::fill_n(blockInfo.begin(), block*2, 2.0);
                    std::fill_n(blockInfo.begin(), block, 1.0);
                    std::fill_n(blockInfo.begin()+(blockInfo.size()-block*2), block*2, 2.0);
                    std::fill_n(blockInfo.begin()+(blockInfo.size()-block), block, 1.0);
                }
                size_t blockIter =0;
                for(size_t j = 0; j < m_chrCount[chr]; ++j){
                    if(m_inclusion[j+index]!=-1){
                        (*snpList).at(m_inclusion[j+index])->SetblockInfo(blockInfo[blockIter]);
                        blockIter++;
                    }
                }

                index+= m_chrCount[chr];
            }
            else{
                throw "Unexpected error";
            }
        }
        else{
            //There is nothing to do with the current chromosome;
            if(m_chrCount.find(chr) != m_chrCount.end()){
                index+= m_chrCount[chr];
            }
        }
    }

}

GenotypeFileHandler::~GenotypeFileHandler(){}


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

			if (bfile_SNP_major) std::cerr << "Detected that binary PED file is v1.00 SNP-major mode" << std::endl;
			else std::cerr << "Detected that binary PED file is v1.00 individual-major mode" << std::endl;

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

void GenotypeFileHandler::skipSnps(size_t const skipNum){
    /** This function is used to skip the chromosome */
    size_t processed = 0;
    while (m_snpIter < m_inputSnp && processed < skipNum){
        size_t indx = 0;
		while ( indx < m_ldSampleSize ){
			std::bitset<8> b; //Initiate the bit array
			char ch[1];
			m_bedFile.read(ch,1); //Read the information
			if (!m_bedFile){
				throw "Problem with the BED file...has the FAM/BIM file been changed?";
			}
			int c = 0;
			while (c<7 && indx < m_ldSampleSize ){
                ++indx;
                c+=2;
			}
		}
		m_snpIter++;
		processed++;
	}
}

ProcessCode GenotypeFileHandler::getSnps(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc, std::vector<Snp*> *snpList, bool &chromosomeStart, bool &chromosomeEnd, double const maf, size_t &prevResidual, size_t &blockSize){
    /** Now that we know exactly how many SNPs are in each chromosome
     *  and how many SNPs that we need from each chromosome, we can
     *  extract the SNPs in a more efficient and easy to understand
     *  manner. Will need to re write this function yet again
     */
    if(m_chrExists.empty()){
        return completed;
    }
    while(m_chrProcessCount.find(m_chrExists.front())==m_chrProcessCount.end()){
        //While we have nothing to do, skip chromosome
        skipSnps(m_chrCount[m_chrExists.front()]);
        //std::cerr << "No SNPs to process, skip chromosome " << m_chrExists.front() << std::endl;
        m_chrExists.pop_front();
        chromosomeStart = true;
        chromosomeEnd = false;
        m_targetProcessed = 0;
        m_processed = 0;
        if(m_chrExists.empty()){
            return completed;
        }
    }
    blockSize = m_blockSizeTract.at(m_chrExists.front());
    size_t processSize = blockSize/3*m_thread;
    if(chromosomeStart) processSize += blockSize/3*2;
    if(m_chrProcessCount[m_chrExists.front()]-(processSize+m_targetProcessed) < blockSize/3){
        //So the last block should be extended
        processSize+=m_chrProcessCount[m_chrExists.front()]-(processSize+m_targetProcessed);
        chromosomeEnd =true;
    }


    prevResidual = genotype.size();
    while(m_snpIter < m_inputSnp){
        bool snp = false;
		if(m_inclusion[m_snpIter] != -1){//indicate whether if we need this snp
			snp=true;
		}
		if(snp){
            genotype.push_back(new Genotype());
            snpLoc.push_back(m_inclusion[m_snpIter]);
		}
		size_t indx = 0; //The iterative count
        double oldM=0.0, newM=0.0,oldS=0.0, newS=0.0;
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
					genotype.back()->AddsampleGenotype(first+second, indx-1); //0 1 or 2
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
			indx > 0 ? genotype.back()->Setmean(newM) : genotype.back()->Setmean(0.0);
			indx > 1 ? genotype.back()->SetstandardDeviation(std::sqrt(newS/(indx - 1.0))) : genotype.back()->SetstandardDeviation(0.0);
            processSize--;
            m_targetProcessed++;
		}
		m_snpIter++;
		m_processed++;
        if(m_processed > m_chrCount[m_chrExists.front()]){
            //Finished this chromosome
            m_chrExists.pop_front(); //Remove the front
			m_processed = 0;
			m_targetProcessed=0;
			//There will not be any case (in my little brain) where we are trying to iterate an empty chromosome because it was dealt with
			chromosomeEnd = true;
            return continueProcess;
        }
        else if(processSize == 0){
            if(m_targetProcessed==m_chrProcessCount[m_chrExists.front()]){
                //Remove the remaining snps
                size_t remaining = m_chrCount[m_chrExists.front()]-m_chrProcessCount[m_chrExists.front()];
                skipSnps(remaining);
                m_chrExists.pop_front(); //Remove the front
                m_processed = 0;
                m_targetProcessed=0;
                chromosomeEnd = true;
            }
            return continueProcess;
        }
    }
    //Now we start processing because there are things to do
    chromosomeEnd=true;
	return completed;

}



ProcessCode GenotypeFileHandler::getSnps(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc, std::vector<Snp*> *snpList, bool &chromosomeStart, bool &chromosomeEnd, double const maf, size_t &numSnp){
	//We will get snps according to the distance
	//We want to use flanking distance, e.g. getting the 1mb flanking on the both side

    size_t processSize = numSnp; //This is the expected number of snps to be processed
	while (m_snpIter < m_inputSnp){ //While there are still Snps to read
		bool snp = false;
		if(m_inclusion[m_snpIter] != -1){//indicate whether if we need this snp
			snp=true;
		}
		if(snp){
            genotype.push_back(new Genotype());
            snpLoc.push_back(m_inclusion[m_snpIter]);
		}
		size_t indx = 0; //The iterative count
        double oldM=0.0, newM=0.0,oldS=0.0, newS=0.0;
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
					if(first == 1 && second == 0) first = 3; //We consider the missing value to 3
					genotype.back()->AddsampleGenotype(first+second, indx-1); //0 1 or 2
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
			indx > 0 ? genotype.back()->Setmean(newM) : genotype.back()->Setmean(0.0);
			indx > 1 ? genotype.back()->SetstandardDeviation(std::sqrt(newS/(indx - 1.0))) : genotype.back()->SetstandardDeviation(0.0);

			double currentMaf = (alleleCount+0.0)/(2*m_ldSampleSize*1.0);
			currentMaf = (currentMaf > 0.5)? 1-currentMaf : currentMaf;
			//remove snps with maf too low or that has 0 variance
			if(maf >= 0.0 && maf > currentMaf){
				std::cerr << "Snp: " << (*snpList)[snpLoc.back()]->GetrsId() << " not included due to maf filtering" << std::endl;
				Genotype *temp = genotype.back();
				genotype.pop_back();
                (*snpList)[snpLoc.back()]->setFlag(0, false);
				snpLoc.pop_back();
				delete temp;

			}
			else{
				processSize--;
			}
		}
		m_snpIter++;
		m_processed++;
		//Check if we have finished the chromosome
        if(m_processed > m_chrCount[m_chrExists.front()]){
			//finished one chromosome
			chromosomeEnd = true;
			m_processed = 0;
			m_chrExists.pop_front();
			return continueProcess;
        }
        //check if we have used all the snps;
        if(processSize == 0){
			return continueProcess;
        }

	}
    chromosomeEnd = true;
	return completed;

}
