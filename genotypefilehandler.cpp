#include "genotypefilehandler.h"

GenotypeFileHandler::GenotypeFileHandler(std::string genotypeFilePrefix, SnpIndex *snpIndex, std::vector<Snp*> &snpList):m_genotypeFilePrefix(genotypeFilePrefix){
	m_chrCount =new SnpIndex();
	m_ldSampleSize = 0;
	m_expectedNumberOfSnp = 0;
	std::vector<std::string> chrCheck; //To check that all the snps are sorted by chromosome
    std::string famFileName = genotypeFilePrefix +".fam";
    std::ifstream famFile(famFileName.c_str(), std::ios::in);
    if(!famFile.is_open()){
        std::cerr << "Cannot open fam file: " << famFileName << std::endl;
        exit(-1);
    }
    std::string line;
    while(std::getline(famFile, line)){
        line = usefulTools::trim(line);
        if(!line.empty()) m_ldSampleSize++;
    }
    famFile.close();
    std::cerr << "A total of " << m_ldSampleSize << " were found in the genotype file for LD construction" << std::endl;
    std::string bimFileName = genotypeFilePrefix+".bim";
    std::ifstream bimFile(bimFileName.c_str(), std::ios::in);
    if(!bimFile.is_open()){
        std::cerr << "Cannot open bim file: " << bimFileName << std::endl;
        exit(-1);
    }
	std::map<std::string, bool> duplicateCheck;
    int duplicateCount = 0;
    while(std::getline(bimFile, line)){
		line = usefulTools::trim(line);
		if(!line.empty()){
            m_expectedNumberOfSnp++;
			std::vector<std::string> token;
            usefulTools::tokenizer(line, "\t ", &token);
            if(token.size() >= 6){
				std::string chr = token[0];
                std::string rs = token[1];
                if(m_chrExists.empty()) m_chrExists.push_back(chr);
                else if(chr.compare(m_chrExists[m_chrExists.size()-1])!= 0) m_chrExists.push_back(chr);
                int snpLoc =-1;
                if(snpIndex->find(rs)){
					snpLoc = snpIndex->value(rs);
					snpList[snpLoc]->setFlag(0, true);
                    if(duplicateCheck.find(rs)!=duplicateCheck.end()){
                        duplicateCount++;
                        m_inclusion.push_back(-1);
                    }
                    else{
                        m_inclusion.push_back(snpLoc);
                        duplicateCheck[rs] = true;
                    }
                }
				if(chrCheck.empty() || chrCheck.size() == 0){
					chrCheck.push_back(chr);
				}
				else if(chrCheck[chrCheck.size()-1].compare(chr)!=0) chrCheck.push_back(chr);
				m_chrCount->increment(chr);
                m_inclusion.push_back(snpLoc);
            }
		}
    }
    bimFile.close();
}

GenotypeFileHandler::~GenotypeFileHandler()
{
	//dtor
}


bool GenotypeFileHandler::openPlinkBinaryFile(const std::string s, std::ifstream & BIT){
	BIT.open(s.c_str(), std::ios::in | std::ios::binary);
	if(!BIT.is_open()){
		std::cerr << "Cannot open the bed file: " << s << std::endl;
		exit(-1);
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
	std::cerr << "check magic number" << std::endl;
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
		std::cerr << " *** Unles you are _sure_ this really is an old BED file *** " << std::endl;
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

