#include "snp.h"

Snp::Snp(std::string chr, std::string rs, size_t bp, size_t sampleSize, std::vector<std::string> &original, std::string refAllele, std::string altAllele, int direction):m_chr(chr), m_rs(rs), m_ref(refAllele), m_alt(altAllele), m_bp(bp), m_sampleSize(sampleSize),m_direction(direction){
    for(size_t i = 0; i < original.size(); ++i){
        if(!usefulTools::isNumeric(original[i])){
            m_original.push_back(0);
            m_remove.push_back(true);
        }
        else if(!std::isfinite(atof(original[i].c_str()))){
            m_original.push_back(0);
            m_remove.push_back(false);
        }
        else{
            m_original.push_back(atof(original[i].c_str()));
            m_remove.push_back(false);
        }
    }
}

void Snp::generateSnpList(boost::ptr_vector<Snp> &snpList, const Command *commander){
    std::ifstream pValue;
    pValue.open(commander->getPvalueFileName().c_str());
    if(!pValue.is_open()){
        throw std::runtime_error("Cannot read the p-value file");
    }
    std::string line;
    //Assume the p-value file should always has a header
    std::getline(pValue, line);
    bool qt = commander->quantitative();
    bool rqt = commander->conRisk();
    bool rcc = commander->diRisk();
    bool isP = commander->isPvalue();
    bool sampleProvided = commander->sampleSizeProvided();


    size_t expectedTokenSize=0;
    size_t chrIndex = commander->getChr();
    expectedTokenSize =(expectedTokenSize<chrIndex)?chrIndex:expectedTokenSize;
    size_t rsIndex = commander->getRs();
    expectedTokenSize =(expectedTokenSize<rsIndex)?rsIndex:expectedTokenSize;
    size_t bpIndex = commander->getBp();
    expectedTokenSize =(expectedTokenSize<bpIndex)?bpIndex:expectedTokenSize;

    size_t sampleIndex = commander->getSampleIndex();
    if(!sampleProvided && (qt || rqt))   expectedTokenSize =(expectedTokenSize<sampleIndex)?sampleIndex:expectedTokenSize;

    //Risk Specific
    size_t dirIndex = commander->getDir();
    if(isP && (rcc || rqt))   expectedTokenSize =(expectedTokenSize<dirIndex)?dirIndex:expectedTokenSize;
    size_t refIndex = commander->getRef();
    if(rcc || rqt)   expectedTokenSize =(expectedTokenSize<refIndex)?refIndex:expectedTokenSize;
    size_t altIndex = commander->getAlt();
    if(rcc || rqt)   expectedTokenSize =(expectedTokenSize<altIndex)?altIndex:expectedTokenSize;
    size_t largestStatIndex =commander->maxStatIndex();
    expectedTokenSize = (expectedTokenSize < largestStatIndex) ?largestStatIndex : expectedTokenSize;
    size_t numIndex = commander->getStatSize();
    size_t samplesize = commander->getSampleSize();

    size_t lineSkipped = 0;

    std::map<std::string, bool> duplication;
    size_t duplicateCount=0;

    while(std::getline(pValue, line)){
        line =usefulTools::trim(line);
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, " \t", &token);
            if(token.size() > expectedTokenSize){
                std::string rsId = token[rsIndex];
                if(duplication.find(rsId) !=duplication.end()){
                    duplicateCount++;
                }
                else{
                    std::string chr = token[chrIndex];
                    size_t bp = atoi(token[bpIndex].c_str());
                    size_t sizeOfSample = samplesize;
                    if(!sampleProvided) sizeOfSample = atoi(token[sampleIndex].c_str());
                    std::string refAllele = "";
                    std::string altAllele = "";
                    int direction = 1;
                    if(rcc || rqt){
                        refAllele = token[refIndex];
                        altAllele = token[altIndex];
                        direction =atof(token[dirIndex].c_str()); //Assumption here, we assume their input is correct
                        if(rcc){
                            //odd ratio
                            direction = (1 <= direction) ?1 : -1;
                        }
                        else if(rqt){
                            //test statistic
                            direction=usefulTools::signum(direction);
                        }
                    }
                    //Now get the statistics
                    std::vector<std::string> statistics;
                    for(size_t i = 0; i < numIndex;++i){
                        statistics.push_back(token[commander->getStatIndex(i)]);
                    }
                    snpList.push_back(new Snp(chr, rsId, bp, sizeOfSample, statistics, refAllele, altAllele,direction));
                    duplication[rsId] = true;
                }
            }
            else{
                lineSkipped++;
            }
        }
    }
    if(lineSkipped!=0){
        std::cerr << lineSkipped << " line(s) skipped as they don't contain enough number of fields." << std::endl;
    }
    pValue.close();
    snpList.sort(Snp::sortSnp);
    std::cerr <<  duplicateCount << " duplicated rsID(s) in the p-value file" << std::endl;
    std::cerr << snpList.size() << " Snps remains" << std::endl;
    if(numIndex != 1) std::cerr << numIndex << " studies processed together" << std::endl;
    if(snpList.size() ==0) throw std::runtime_error("Programme terminated as there are no snp provided");
}


bool Snp::sortSnp (const Snp& i, const Snp& j){
    if(i.getChr().compare(j.getChr()) == 0)
		if(i.getBp() == j.getBp())
			return i.getRs().compare(j.getRs()) < 0;
		else
			return i.getBp() < j.getBp();
	else return (i.getChr().compare(j.getChr()) < 0);
}




