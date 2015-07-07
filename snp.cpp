#include "snp.h"

size_t Snp::m_maxSampleSize=0;
size_t Snp::m_perfectId=0;
double Snp::m_adjustment = 1.0;

Snp::Snp(std::string chr, std::string rs, size_t bp, double sampleSize, double original):m_chr(chr), m_rs(rs), m_bp(bp), m_sampleSize(sampleSize), m_original(original){
	m_beta = std::make_shared<double>(original);
	m_sqrtChiSq = std::make_shared<double>(original);
	m_heritability = std::make_shared<double>(0.0);
	m_effectiveNumber=0.0;
	m_additionVariance =0.0;
	m_variance =0.0;
	m_targetClass = this;
	m_sign = usefulTools::signum(original);
	m_perfectLdId= 0;
	if(Snp::m_maxSampleSize < m_sampleSize){
        Snp::m_maxSampleSize = m_sampleSize;
	}
}

std::string Snp::Getchr() const { return m_chr; }
std::string Snp::GetrsId() const { return m_rs; }
size_t Snp::GetperfectId() const { return m_perfectLdId; }
size_t Snp::Getbp() const { return m_bp; }
size_t Snp::GetregionSize() const {return m_regionFlag.size(); }
size_t Snp::GetmaxSampleSize() {return Snp::m_maxSampleSize; }
double Snp::GetsampleSize() const { return m_sampleSize; }
double Snp::Getoriginal() const { return m_original; }
double Snp::Getadjustment() {return Snp::m_adjustment; }
double Snp::Getbeta() const {
	return (*m_beta)/(double)(m_beta.use_count());
}
double Snp::GetsignedSqrtChiSq() const {
	return (*m_sqrtChiSq)/(double)(m_sqrtChiSq.use_count());
}

void Snp::Setheritability(double heritability ) { (*m_heritability) = heritability;}
void Snp::Setvariance(double i){ m_variance = i; }
void Snp::SetadditionVariance(double i){ m_additionVariance = i; }
void Snp::Setvariance(double const sigma, double const sigmaSquared, double const sigmaPowerThree, double const sigmaPowerFour ){
    m_variance = sigma;
    m_additionVariance = sigmaSquared;
    m_sigmaPowerThree = sigmaPowerThree;
    m_sigmaPowerFour = sigmaPowerFour;
}

void Snp::Setsign(int directionEffect){ m_sign = directionEffect; }
double Snp::GetheritabilityChi() const { return (*m_heritability)/(double)(m_beta.use_count()); }
double Snp::Getvariance() const { return m_variance; }
double Snp::Getheritability() const {
	double heritability = (*m_heritability)/(double)(m_beta.use_count());
    heritability = heritability*heritability-1.0/(m_sampleSize*1.0);
    //std::cout << m_rs << "\tSample size: " << m_sampleSize << std::endl;
	return heritability;
}

void Snp::Setadjustment(const double prevalence, const size_t caseSize, const size_t controlSize){
    double portionCase = (caseSize+0.0) / (caseSize+controlSize+0.0);
	double i2 = usefulTools::dnorm(usefulTools::qnorm(prevalence))/(prevalence);
	i2 = i2*i2;
    Snp::m_adjustment = ((1-prevalence)*(1-prevalence))/(i2*portionCase*(1-portionCase));
}


void Snp::shareHeritability( Snp* i ){
	if(i->m_beta == m_beta){
		return;
	}
	if(i->m_perfectLdId != 0) m_perfectLdId = i->m_perfectLdId;
	else if(m_perfectLdId != 0) i->m_perfectLdId = m_perfectLdId;
	else if(m_perfectLdId != 0 && i->m_perfectLdId != 0 && m_perfectLdId != i->m_perfectLdId) throw "Unexpected behaviour, different perfectLd id detected";
    else if(m_perfectLdId != 0 && i->m_perfectLdId != 0){/**do nothing */}
    else if(m_perfectLdId == 0 && i->m_perfectLdId == 0){
        m_perfectLdId = Snp::m_perfectId;
        i->m_perfectLdId = Snp::m_perfectId;
        Snp::m_perfectId++;
    }
	(*i->m_beta) += (*m_beta); //Add up the beta to get an average
	(*i->m_sqrtChiSq) += (*m_sqrtChiSq);
	(*i->m_heritability) += (*m_heritability);
	//i is the one who buy
	//this is the one who sell
    //For this's family, all of them should point to the same location
	m_beta = (i->m_beta); //They now share the same beta
	m_sqrtChiSq = (i->m_sqrtChiSq);
	m_heritability = (i->m_heritability);
	Snp* currentClass = m_targetClass;
	Snp* prevClass = this;
	while(currentClass != this){
		currentClass->m_beta = i->m_beta; //The subsequent stuff are also pointing here
		currentClass->m_sqrtChiSq = i->m_sqrtChiSq;
		currentClass->m_heritability = i->m_heritability;
        prevClass = currentClass;
        currentClass = currentClass->m_targetClass;
    }
    //Now they are all pointing to the same value
	//Now prevClass is the last person in chain
    prevClass->m_targetClass = i->m_targetClass;
    i->m_targetClass = this;

}


bool Snp::GetFlag(size_t index) const {
	if(index>= m_regionFlag.size()) return false;
	return m_regionFlag[index];
}


Snp::~Snp(){}

void Snp::generateSnpList(std::vector<Snp*> &snpList, const Command *commander){
	std::ifstream pValue;
    pValue.open(commander->GetpValueFileName().c_str());
    if(!pValue.is_open()){
        throw "Cannot read the p-value file";
    }
    std::string line;
    if(commander->hasHeader()){ //If the p-value file contain a header, we will ignore it
        std::getline(pValue, line);
    }
    size_t bpIndex = commander->GetbpIndex();
    size_t chrIndex = commander->GetchrIndex();
    size_t rsIndex = commander->GetrsIndex();
    size_t index = 0, sIndex=0;
    if(commander->caseControl()) index=commander->GetcIndex();
    if(commander->quantitative()) index=commander->GettIndex();
    if(!commander->provideSampleSize()) sIndex= commander->GetsampleSizeIndex();
    while(std::getline(pValue, line)){
        line =usefulTools::trim(line);
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, " \t", &token);
            if(token.size() > bpIndex && token.size() > chrIndex &&
               token.size() > rsIndex && token.size() > index &&
               token.size() > sIndex){
                std::string chr = token[chrIndex];
                size_t bp = atoi(token[bpIndex].c_str());
                std::string rsId = token[rsIndex];
                size_t sizeOfSample = commander->GetsampleSize();
                if(!commander->provideSampleSize()) sizeOfSample = atoi(token[sIndex].c_str());
                double predictedBeta = atof(token[index].c_str());
                if(!std::isfinite(predictedBeta)){
                    std::cerr << rsId << " does not have finite input ("<< predictedBeta <<"), will set it to zero" << std::endl;
                    predictedBeta=0;
                }
                snpList.push_back(new Snp(chr, rsId, bp, sizeOfSample, predictedBeta));
            }
        }
    }
    pValue.close();
    std::sort(snpList.begin(), snpList.end(), Snp::sortSnp);
    snpList.erase( unique( snpList.begin(), snpList.end() ), snpList.end() );
    std::cerr << "There are a total of " << snpList.size() << " Snps in the input" << std::endl;
    if(snpList.size() ==0) throw "Programme terminated as there are no snp provided";
}

bool Snp::sortSnp (Snp* i, Snp* j){
    if(i->Getchr().compare(j->Getchr()) == 0)
		if(i->Getbp() == j->Getbp())
			return i->GetrsId().compare(j->GetrsId()) < 0;
		else
			return i->Getbp() < j->Getbp();
	else return (i->Getchr().compare(j->Getchr()) < 0);
}


void Snp::generateSnpIndex(SnpIndex *snpIndex, std::vector<Snp*> &snpList, Region *regionList, bool isPvalue, double extremeRatio){
	std::vector<size_t> regionIncrementationIndex(regionList->GetnumRegion(), 0);
	size_t duplicate = 0;
	for(size_t i = 0; i < snpList.size(); ++i){
        //If the snp is new
        if(!snpIndex->contains(snpList[i]->GetrsId())){
            snpIndex->set(snpList[i]->GetrsId(), i);
            try{
                snpList[i]->computeVarianceExplainedChi(isPvalue, extremeRatio);
                //The default flag (with LD), is always false at this stage
                snpList[i]->m_regionFlag.push_back(false);

                if(regionList->GetnumRegion() != 0){
                    std::vector<bool> padding(regionList->GetnumRegion(), false);
                    snpList[i]->m_regionFlag.insert(snpList[i]->m_regionFlag.end(), padding.begin(), padding.end());
                }
                for(size_t j = 0; j < regionList->GetnumRegion(); ++j){
                    for(unsigned k = regionIncrementationIndex[j]; k < regionList->GetintervalSize(j); ++k){
                        //check whether if this snp falls within the region
                        if(regionList->Getchr(j,k).compare(snpList[i]->Getchr())==0 &&
                           regionList->Getstart(j,k) <= snpList[i]->Getbp() &&
                           regionList->Getend(j,k) >= snpList[i]->Getbp()){
                            regionIncrementationIndex[j] = k;
                            snpList[i]->setFlag(j+1, true);
                            break;
                           }
                    }
                }
            }
            catch(const char* e){
                std::cerr << e << std::endl;
            }

        }
        else{
            duplicate++;
        }

    }
    if(duplicate == 0) std::cerr << "There are no duplicated rsID in the p-value file" << std::endl << std::endl;
    else std::cerr <<  "There are a total of " << duplicate << " duplicated rsID(s) in the p-value file" << std::endl << std::endl;
}

void Snp::addDirection(SnpIndex *snpIndex, std::vector<Snp*> &snpList,std::string dirFile){
    if(dirFile.empty()){
        throw "Should not happen as we doesn't require the direction file";
    }
    std::ifstream direction;
    direction.open(dirFile.c_str());
    if(!direction.is_open()){
        throw "Cannot open the direction file";
    }
    //Assume no header
    std::string line;
    while(std::getline(direction, line)){
        line = usefulTools::trim(line);
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, "\t ", &token);
            //The first should be rsId, the second should be the direction
            if(token.size() > 1){
                if(snpIndex->contains(token[0])){
                    size_t refId = snpIndex->value(token[0]);
                    snpList[refId]->Setsign(atoi(line.c_str()));
                }
                else{
                    throw "The direction file contains Snps not found in p-value file. Most likely they are not matched. Please check your input!";
                }
            }
        }
    }
    direction.close();
}

void Snp::generateSnpIndex(SnpIndex *snpIndex, std::vector<Snp*> &snpList, const size_t &caseSize, const size_t &controlSize, const double &prevalence, Region *regionList, bool isPvalue){
	std::vector<size_t> regionIncrementationIndex(regionList->GetnumRegion(), 0);
	size_t duplicate = 0;
	for(size_t i = 0; i < snpList.size(); ++i){
        if(!snpIndex->contains(snpList[i]->GetrsId())){
            snpIndex->set(snpList[i]->GetrsId(), i);
            try{
                snpList[i]->computeVarianceExplainedChi(caseSize, controlSize, prevalence, isPvalue);;
                snpList[i]->m_regionFlag.push_back(false);
                if(regionList->GetnumRegion() != 0){
                    std::vector<bool> padding(regionList->GetnumRegion(), false);
                    snpList[i]->m_regionFlag.insert(snpList[i]->m_regionFlag.end(), padding.begin(), padding.end());
                }
                for(size_t j = 0; j < regionList->GetnumRegion(); ++j){
                    for(unsigned k = regionIncrementationIndex[j]; k < regionList->GetintervalSize(j); ++k){
                        if(regionList->Getchr(j,k).compare(snpList[i]->Getchr())==0 &&
                           regionList->Getstart(j,k) <= snpList[i]->Getbp() &&
                           regionList->Getend(j,k) >= snpList[i]->Getbp()){
                            regionIncrementationIndex[j] = k;
                            snpList[i]->setFlag(j+1, true);
                            break;
                           }
                    }
                }
            }
            catch(const char * e){
                std::cerr << e << std::endl;
            }
        }
        else{
            duplicate++;
        }

    }
    if(duplicate == 0) std::cerr << "There are no duplicated rsID in the p-value file" << std::endl << std::endl;
    else std::cerr <<  "There are a total of " << duplicate << " duplicated rsID(s) in the p-value file" << std::endl << std::endl;
}

void Snp::computeVarianceExplainedChi(bool isPvalue, double extremeRatio){
    if(isPvalue){
        (*m_beta) = usefulTools::qnorm(1.0-((m_original+0.0)/2.0));
        if(!std::isfinite((*m_beta))) (*m_beta) = usefulTools::qnorm(((m_original+0.0)/2.0));
        if(!std::isfinite((*m_beta))){
            //This will only happen when the p-value is either 1 or 0
            if(m_original >= 1.0){
                (*m_beta) = 0.0; //There is no effect anyway
            }
            else{
                throw "WARNING! A p-value of 0 is observed. We don't know how to convert it into chi square. Will skip this snp "+m_rs;
            }
        }
    }
    (*m_beta) = (*m_beta)*(*m_beta);
    (*m_sqrtChiSq)= sqrt((*m_beta))*m_sign;
	(*m_beta) = ((*m_beta)-1)/(m_sampleSize-2.0+(*m_beta));
    Snp::m_adjustment = extremeRatio;

}

void Snp::computeVarianceExplainedChi(const size_t &caseSize, const size_t &controlSize, const double &prevalence, bool isPvalue){
    if(isPvalue){
        (*m_beta) = usefulTools::qnorm(1.0-((m_original+0.0)/2.0));
        if(!std::isfinite((*m_beta))) (*m_beta) =usefulTools::qnorm(((m_original+0.0)/2.0));
        if(!std::isfinite((*m_beta))){
            //This will only happen when the p-value is either 1 or 0
            if(m_original >= 1.0){
                (*m_beta) = 0.0; //There is no effect anyway
            }
            else{
                throw "WARNING! A p-value of 0 is observed. We don't know how to convert it into chi square. Will skip this snp "+m_rs;
            }
        }
        (*m_beta) = (*m_beta)*(*m_beta);
    }
    (*m_sqrtChiSq)= sqrt((*m_beta))*m_sign;
	int totalSampleSize = caseSize + controlSize;
    (*m_beta) = ((*m_beta)-1.0)/(totalSampleSize -2.0+(*m_beta));
    Snp::m_maxSampleSize = totalSampleSize;
}

void Snp::setFlag(size_t index, bool value){
	if(index >= m_regionFlag.size()){
		throw "Error: region flag out of bound!";
	}
    m_regionFlag[index] = value;
}


void Snp::cleanSnp(std::vector<Snp*> &snpList){
    for(size_t i= 0; i < snpList.size(); ++i){
        delete snpList[i];
    }
    snpList.clear();
}

bool Snp::Concordant(std::string chr, size_t bp, std::string rsId) const{
    return chr.compare(m_chr) ==0 && bp==m_bp && rsId.compare(m_rs) == 0;
}






