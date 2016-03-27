#include "genotypefilehandler.h"

GenotypeFileHandler::GenotypeFileHandler()
{
    //ctor
}

GenotypeFileHandler::~GenotypeFileHandler()
{
    //dtor
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
    m_genotypeFilePrefix = commander.getReferenceFilePrefix();
    m_mafThreshold = commander.getMAF();
    m_blockSize = commander.getSizeOfBlock();

    //First check if the bed file format is correct, if not, we just stop
    std::string bedFileName = m_genotypeFilePrefix+".bed";
    bool bfile_SNP_major = openPlinkBinaryFile(bedFileName, m_bedFile); //We will try to open the connection to bedFile
    if(!bfile_SNP_major) throw std::runtime_error("We currently have no plan of implementing the individual-major mode. Please use the snp-major format");

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

    std::string bimFileName = m_genotypeFilePrefix+".bim";
    std::ifstream bimFile;
    bimFile.open(bimFileName.c_str());
    if(!bimFile.is_open())  throw std::runtime_error("Cannot open bim file");

    size_t emptyLineCheck = 0; // Only use for checking for empty line, consider this as a safeguard
    std::map<std::string, bool> duplicateCheck;
    size_t nDuplicated=0,nInvalid=0, nAmbig=0;
    while(std::getline(bimFile, line)){
        line = usefulTools::trim(line);
        //The only thing I cannot filter beforehand is the MAF
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, "\t", &token);
//            m_inclusion.push_back(-1);
            if(token.size()>=6){
                std::string rs = token[1];
                if(duplicateCheck.find(rs)!=duplicateCheck.end())   nDuplicated++; // Don't want it if it is duplicated (impossible here actually)
                else if(snpIndex.find(rs)!=snpIndex.end()){
                    size_t location = snpIndex.at(rs);
                    bool ambig = false;
                    if(snpList.at(location).concordant(token[0], atoi(token[3].c_str()), rs, token[4], token[5], ambig )){
                        if(!m_keepAmbiguous && ambig) nAmbig ++; // We don't want to keep ambiguous SNPs and this is an ambiguous SNP
                        else{
                            m_snpLineNumber.push_back(m_nSnp);
//                            m_inclusion.back()=location;
                            m_inclusion.push_back(location);
                            snpList[location].setFlag(0,true);
                        }
                    }
                    else nInvalid++;
                }
            }
            else throw "Malformed bim file, does not contain all necessary information";
            m_nSnp++;
        }
        else emptyLineCheck++;
    }
    m_nRequiredSnps = m_inclusion.size();
    m_nBytes=ceil((double)m_nRefSample/4.0);
    bimFile.close();
    fprintf(stderr,"\nReference Panel Information:\n");
    fprintf(stderr,"==================================\n");
    fprintf(stderr, "%lu samples were found in the reference panel\n", m_nRefSample);
    fprintf(stderr, "%lu discordant SNP(s) identified\n", nInvalid);
    if(nAmbig!=0)
        fprintf(stderr, "%lu ambiguous SNP(s) filtered\n", nAmbig);
    if(nDuplicated)
        fprintf(stderr, "%lu duplicated SNP(s) found in reference panel\n", nDuplicated);
    fprintf(stderr, "%lu SNP(s) will be included\n", m_nRequiredSnps);
    m_snpIter=0;
}



// Find the first SNP that for us to include.
// Then read all the SNPs within the region

void GenotypeFileHandler::getBlock(boost::ptr_vector<Snp> &snpList, boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, bool &windowEnd, bool &completed, std::vector<size_t> &boundary, bool addition){
    if(m_snpIter >= m_nRequiredSnps) return; //Nothing to read
    bool starting  = !addition;// Indicate that we are trying to identify the start of the current block
    // Simple declaration
    int blockStartIndex = (addition)? snpLoc[boundary.back()] : 0; // This is index for the start SNP of this block
    std::string blockChr= (addition)? snpList[blockStartIndex].getChr(): ""; // Checking if there is any chromosome change
    size_t blockStartLoc= (addition)? snpList[blockStartIndex].getLoc(): 0; // The bp of the first SNP included in this block
    size_t lastUsedLoc = (addition)? snpList[snpLoc.back()].getLoc() : 0; // The bp of the last SNP included in this block
    //m_inclusion contains the snpLoc whereas the m_snpLineNumber tells us how many line do we have to skip
    for(; m_snpIter < m_nRequiredSnps; ++m_snpIter){
        int snpIndex = m_inclusion[m_snpIter];
        if(starting){
            blockChr = snpList[snpIndex].getChr();
            blockStartLoc = snpList[snpIndex].getLoc();
            lastUsedLoc = blockStartLoc;
        }
        size_t currentLoc = snpList[snpIndex].getLoc();
        std::string currentChr = snpList[snpIndex].getChr();
        if(currentChr.compare(blockChr)!=0 ||
            (currentLoc-lastUsedLoc) > m_blockSize)
        {
                windowEnd=true;
                return;
        }
        else if(currentLoc-blockStartLoc > m_blockSize) return;
        else{
            Genotype *tempGenotype = new Genotype();
            size_t indx = 0;
            size_t alleleCount=0, validSample = 0;
            m_bedFile.seekg((m_snpLineNumber[m_snpIter]-m_nSnpSkipped) * m_nBytes, m_bedFile.cur); //Now read the target genotype
            m_nSnpSkipped = m_snpLineNumber[m_snpIter]+1;
            std::bitset<8> byteHolder; //Initiate the bit array
            char allGeno[m_nBytes]; //This reads in all the samples into ch
            m_bedFile.read(allGeno, m_nBytes);
            if (!m_bedFile) throw std::runtime_error("Problem with the BED file...has the FAM/BIM file been changed?");
            size_t sampleIndex = 0;
            for(size_t byteRunner=0; byteRunner < m_nBytes; ++byteRunner){ //Each byte contain 4 samples, thus m_nBytes = sampleSize/4
                byteHolder=allGeno[byteRunner];
                size_t genoBit = 0;
                while (genoBit<7 && sampleIndex<m_nRefSample ){
                    int first = byteHolder[genoBit++];
                    int second = byteHolder[genoBit++];
                    if(!(first==1 && second==0)){ //When not missing
                        validSample++;
                        alleleCount += first+second;
                    }
                    tempGenotype->AddsampleGenotype(first,second, sampleIndex);
                    sampleIndex++;
                }
            }
            double currentMaf = (validSample>0)?((double)alleleCount)/(2.0*(double)validSample):-1.0;
            currentMaf = (currentMaf > 0.5)? 1-currentMaf : currentMaf;
            if(m_mafThreshold > currentMaf){
                m_inclusion[m_snpIter]=-1;
                m_nFilter++;
                snpList[snpIndex].setFlag(0,false);
                snpList[snpIndex].setStatus('F');
                delete tempGenotype;
            }
            else{
                if(starting){
                    boundary.push_back(snpLoc.size()); // If this is the start, put in the snpLoc index
                    starting = false;
                }
                snpList[snpIndex].setStatus('I');
                genotype.push_back(tempGenotype);
                snpLoc.push_back(snpIndex); // index of in snpList
                lastUsedLoc  = currentLoc; // This is for the coordinates
            }
        }
    }
    windowEnd = true;
    completed = true;
    m_bedFile.close();
}
