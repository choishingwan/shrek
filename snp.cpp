#include "snp.h"
/*
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
*/
Snp::Snp(std::string chr, std::string rs, size_t bp, size_t sampleSize, double original, std::string refAllele, std::string altAllele, int direction):m_chr(chr), m_rs(rs), m_ref(refAllele), m_alt(altAllele), m_bp(bp), m_sampleSize(sampleSize),m_original(original),m_direction(direction){}
Snp::~Snp(){}

void Snp::computeVarianceExplained(const Command *commander){
    //There are 4 possibilities
    bool qt = commander->quantitative();
    bool cc = commander->caseControl();
    bool rqt = commander->conRisk();
    bool rcc = commander->diRisk();
    bool isP = commander->isPvalue();
    if(qt || rqt){
        //Sample size should already be obtained for the qt runs
        if(isP){
            double beta = usefulTools::qnorm(1.0-((m_original+0.0)/2.0));
            if(m_original >= 1.0) beta = 0.0;
            else if(m_original ==0.0){ //This is unsafe as == of double is always a problem
                beta=0.0;
            }
            beta = beta*m_direction;
            if(rqt) beta = beta/sqrt(m_sampleSize-2.0+beta*beta);
            else{
                beta = beta*beta;
                beta = (beta-1.0)/(m_sampleSize-2.0+beta);
            }
            m_beta=beta;
            m_heritability=0.0;
        }
        else{
            double beta = m_original;
            if(rqt) beta = beta/sqrt(m_sampleSize-2.0+beta*beta);
            else if(qt){
                beta=beta*beta;
                beta = (beta-1.0)/(m_sampleSize-2.0+beta);
            }
            m_beta=beta;
            m_heritability=0.0;
        }
    }
    else if(cc || rcc){
        size_t caseSize=commander->getCaseSize();
        size_t controlSize=commander->getControlSize();
        if(isP){
            double beta = usefulTools::qnorm(1.0-((m_original+0.0)/2.0));
            if(m_original >= 1.0) beta = 0.0;
            else if(m_original ==0.0){//This is unsafe as == of double is always a problem
                beta=0.0;
            }
            beta = beta*m_direction;
            if(rcc) beta = beta/sqrt(caseSize+controlSize-2.0+beta*beta);
            else{
                beta = beta*beta;
                beta = (beta-1.0)/(caseSize+controlSize-2.0+beta);
            }
            m_beta=beta;
            m_heritability=0.0;
        }
        else{
            double beta =m_original;
            if(rcc){
                beta = sqrt(beta)*m_direction;
                beta = (beta)/(caseSize+controlSize -2.0+beta*beta);
            }
            else if(cc){
                beta = (beta-1.0)/(caseSize+controlSize -2.0+beta);
            }
            m_beta=beta;
            m_heritability=0.0;
        }
    }
    else{
        throw std::runtime_error("Undefined mode for SNP processing");
    }
}


void Snp::generateSnpIndex(std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, const Command *commander, Region *regionList){
    std::vector<size_t> regionIncrementationIndex(regionList->getNumRegion(), 0);
	size_t duplicate = 0;

	for(size_t i = 0; i < snpList.size(); ++i){
        //If the snp is new
        if(snpIndex.find(snpList[i].getRs())== snpIndex.end()){
            snpIndex[snpList[i].getRs()] =i ;

            snpList[i].computeVarianceExplained(commander);
            //The default flag (with LD), is always false at this stage
            snpList[i].m_regionFlag.push_back(false);
            if(regionList->getNumRegion() != 0){
                std::vector<bool> padding(regionList->getNumRegion(), false);
                snpList[i].m_regionFlag.insert(snpList[i].m_regionFlag.end(), padding.begin(), padding.end());
            }
            for(size_t j = 0; j < regionList->getNumRegion(); ++j){
                for(unsigned k = regionIncrementationIndex.at(j); k < regionList->getIntervalSize(j); ++k){
                    //check whether if this snp falls within the region
                    if( regionList->getChr(j,k).compare(snpList[i].getChr())==0 &&
                        regionList->getStart(j,k) <= snpList[i].getBp() &&
                        regionList->getEnd(j,k) >= snpList[i].getBp()){

                        regionIncrementationIndex.at(j) = k;
                        snpList[i].setFlag(j+1, true);
                        break;
                    }
                }
            }
        }
        else{
            duplicate++;
        }

    }
}

void Snp::setFlag(const size_t i, bool flag){
    m_regionFlag.at(i) = flag;
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
    //size_t largestStatIndex =commander->maxStatIndex();
    //expectedTokenSize = (expectedTokenSize < largestStatIndex) ?largestStatIndex : expectedTokenSize;
    size_t statIndex = commander->getStat();
    expectedTokenSize = (expectedTokenSize < statIndex) ?statIndex : expectedTokenSize;

    //size_t numIndex = commander->getStatSize();
    size_t samplesize = commander->getSampleSize();

    size_t lineSkipped = 0;

    std::map<std::string, bool> duplication;
    size_t duplicateCount=0;
    size_t removeCount = 0;
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
                    //std::vector<std::string> statistics;
                    //for(size_t i = 0; i < numIndex;++i){
                    //    statistics.push_back(token[commander->getStatIndex(i)]);
                    //}
                    if(!usefulTools::isNumeric(token[statIndex])){
                        removeCount++;
                    }
                    else{
                        double statistic = atof(token[statIndex].c_str());
                        if(statistic == 0.0 && isP){ //This is unsafe as == of double is always a problem
                            removeCount++;
                        }
                        else{
                            snpList.push_back(new Snp(chr, rsId, bp, sizeOfSample, statistic, refAllele, altAllele,direction));
                            duplication[rsId] = true;
                        }
                    }
                }
            }
            else{
                lineSkipped++;
            }
        }
    }
    if(lineSkipped!=0){
        std::cerr << lineSkipped << " line(s) skipped." << std::endl;
    }
    pValue.close();
    snpList.sort(Snp::sortSnp);
    std::cerr <<  duplicateCount << " duplicated rsID(s) in the p-value file" << std::endl;
    std::cerr << removeCount << " Snps with p-value = 0 or non-numeric statistics removed" << std::endl;
    std::cerr << snpList.size() << " Snps remains" << std::endl;
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

