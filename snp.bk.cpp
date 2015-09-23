
#include "snp.h"

size_t Snp::m_maxSampleSize=0;
size_t Snp::m_perfectId=0;
double Snp::m_adjustment = 1.0;

Snp::Snp(std::string chr, std::string rs, size_t bp, double sampleSize, double original, std::string refAllele, std::string altAllele):m_chr(chr), m_rs(rs), m_bp(bp), m_sampleSize(sampleSize), m_original(original){
	m_beta = std::make_shared<double>(original);
	m_sqrtChiSq = std::make_shared<double>(original);
	m_heritability = std::make_shared<double>(0.0);
	m_effectiveNumber=0.0;
	m_snpLDSC = 0.0;
	m_targetClass = this;
	m_sign = usefulTools::signum(original);
	m_perfectLdId= 0;
	if(Snp::m_maxSampleSize < m_sampleSize){
        Snp::m_maxSampleSize = m_sampleSize;
	}
    m_ref = refAllele;
    m_alt = altAllele;
}

std::string Snp::Getchr() const { return m_chr; }
std::string Snp::GetrsId() const { return m_rs; }
std::string Snp::Getref() const { return m_ref; }
std::string Snp::Getalt() const { return m_alt; }
size_t Snp::GetperfectId() const { return m_perfectLdId; }
size_t Snp::Getbp() const { return m_bp; }
size_t Snp::GetregionSize() const {return m_regionFlag.size(); }
size_t Snp::GetblockInfo() const {return m_blockInfo; }
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
double Snp::GeteffectiveNumber() const { return m_effectiveNumber; }
double Snp::GetsnpLDSC() const { return m_snpLDSC; }
void Snp::Setheritability(double heritability ) {
     (*m_heritability) = heritability;
}
void Snp::SeteffectiveNumber(double effective){
    m_effectiveNumber = effective;
}
void Snp::SetsnpLDSC(double effect){
    m_snpLDSC = effect;
}
void Snp::Setvariance(double i){ m_variance = i; }
void Snp::SetblockInfo(size_t blockInfo){m_blockInfo = blockInfo;}
void Snp::SetadditionVariance(double i){ m_additionVariance = i; }
void Snp::Setvariance(double const sigma, double const sigmaSquared, double const sigmaPowerThree, double const sigmaPowerFour ){
    m_variance = sigma;
    m_additionVariance = sigmaSquared;
    m_sigmaPowerThree = sigmaPowerThree;
    m_sigmaPowerFour = sigmaPowerFour;
}

void Snp::Setsign(int directionEffect){ m_sign = directionEffect; }
double Snp::Getheritability() const { return (*m_heritability)/(double)(m_beta.use_count()); }
double Snp::Getvariance() const { return m_variance; }

void Snp::Setadjustment(const double prevalence, const size_t caseSize, const size_t controlSize){
    double portionCase = (double)(caseSize) / (double)(caseSize+controlSize);
	double i2 = usefulTools::dnorm(usefulTools::qnorm(prevalence))/(prevalence);
	i2 = i2*i2;
    Snp::m_adjustment = ((1.0-prevalence)*(1.0-prevalence))/(i2*portionCase*(1-portionCase));
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
	return m_regionFlag.at(index);
}


Snp::~Snp(){}

void Snp::generateSnpList(std::vector<Snp*> &snpList, const Command *commander){
	std::ifstream pValue;
    pValue.open(commander->GetpValueFileName().c_str());
    if(!pValue.is_open()){
        throw "Cannot read the p-value file";
    }
    std::string line;
    //Assume the p-value file should always has a header
    std::getline(pValue, line);

    size_t bpIndex = commander->GetbpIndex();
    size_t maxIndex= bpIndex;
    size_t chrIndex = commander->GetchrIndex();
    if(chrIndex > maxIndex) maxIndex =chrIndex;
    size_t rsIndex = commander->GetrsIndex();
    if(rsIndex > maxIndex) maxIndex = rsIndex;
    size_t sIndex = 0;
    size_t index = commander->GetIndex();
    if(index > maxIndex) maxIndex = index;
    size_t altIndex = 0;
    size_t refIndex = 0;
    std::map<std::string, bool> duplication;
    size_t duplicateCount=0;
    if(commander->risk()){
        altIndex = commander->GetaltIndex();
        if(altIndex > maxIndex) maxIndex = maxIndex;
        refIndex = commander->GetrefIndex();
        if(refIndex > maxIndex) maxIndex = refIndex;
    }
    if(!commander->provideSampleSize()) sIndex= commander->GetsampleSizeIndex();
    if(sIndex > maxIndex) maxIndex = sIndex;
    std::string removeSnps="";
    size_t removeSnpCount = 0;
    while(std::getline(pValue, line)){
        line =usefulTools::trim(line);
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, " \t", &token);
            if(token.size() > maxIndex){
                std::string chr = token[chrIndex];
                size_t bp = atoi(token[bpIndex].c_str());
                std::string rsId = token[rsIndex];
                size_t sizeOfSample = 0;
                if(duplication.find(rsId) !=duplication.end()){
                    duplicateCount++;
                }
                else{
                    if(!commander->provideSampleSize()) sizeOfSample = atoi(token[sIndex].c_str());
                    else sizeOfSample = commander->GetsampleSize();
                    std::string refAllele = "";
                    std::string altAllele = "";
                    if(commander->risk()){
                        refAllele = token[refIndex];
                        altAllele = token[altIndex];

                    }
                    if(!usefulTools::isNumeric(token[index])){
                        //Check if the input is a number. If it is not, then it should be filtered out.
                        removeSnps.append(rsId);
                        removeSnps.append(",");
                        //std::cerr << "Remove: " << rsId << "\t" << token[index] << std::endl;
                        removeSnpCount++;
                    }
                    else{
                        double predictedBeta = atof(token[index].c_str());
                        if(!std::isfinite(predictedBeta)){
                            std::cerr << rsId << " does not have finite input ("<< predictedBeta <<"), will set it to zero" << std::endl;
                            predictedBeta=0;
                        }
                        snpList.push_back(new Snp(chr, rsId, bp, sizeOfSample, predictedBeta, refAllele, altAllele));
                        duplication[rsId] = true;
                    }
                }
            }
        }
    }
    if(!removeSnps.empty()){
        std::cerr << removeSnpCount << " SNPs removed with missing values." << std::endl;
    }
    pValue.close();

    //This can be slow when we have a large amount of SNPs. Remove all duplicated SNPs here might be a good choice
    std::sort(snpList.begin(), snpList.end(), Snp::sortSnp);
//    snpList.erase( unique( snpList.begin(), snpList.end() ), snpList.end() );
    if(duplicateCount == 0) std::cerr << "There are no duplicated rsID in the p-value file" << std::endl << std::endl;
    else std::cerr <<  "There are a total of " << duplicateCount << " duplicated rsID(s) in the p-value file" << std::endl << std::endl;
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


void Snp::generateSnpIndex(std::map<std::string, size_t> &snpIndex, std::vector<Snp*> &snpList, Region *regionList, bool isPvalue, double extremeRatio){
	std::vector<size_t> regionIncrementationIndex(regionList->GetnumRegion(), 0);
	size_t duplicate = 0;
	for(size_t i = 0; i < snpList.size(); ++i){
        //If the snp is new
        if(snpIndex.find(snpList[i]->GetrsId())== snpIndex.end()){
            snpIndex[snpList[i]->GetrsId()] =i ;
            try{
                snpList[i]->computeVarianceExplainedChi(isPvalue, extremeRatio);
                //The default flag (with LD), is always false at this stage
                snpList[i]->m_regionFlag.push_back(false);

                if(regionList->GetnumRegion() != 0){
                    std::vector<bool> padding(regionList->GetnumRegion(), false);
                    snpList[i]->m_regionFlag.insert(snpList[i]->m_regionFlag.end(), padding.begin(), padding.end());
                }
                for(size_t j = 0; j < regionList->GetnumRegion(); ++j){
                    for(unsigned k = regionIncrementationIndex.at(j); k < regionList->GetintervalSize(j); ++k){
                        //check whether if this snp falls within the region
                        if(regionList->Getchr(j,k).compare(snpList[i]->Getchr())==0 &&
                           regionList->Getstart(j,k) <= snpList[i]->Getbp() &&
                           regionList->Getend(j,k) >= snpList[i]->Getbp()){
                            regionIncrementationIndex.at(j) = k;
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

}


void Snp::generateSnpIndex(std::map<std::string, size_t> &snpIndex, std::vector<Snp*> &snpList, const size_t &caseSize, const size_t &controlSize, const double &prevalence, Region *regionList, bool isPvalue){
	std::vector<size_t> regionIncrementationIndex(regionList->GetnumRegion(), 0);
	size_t duplicate = 0;
	for(size_t i = 0; i < snpList.size(); ++i){
        if(snpIndex.find(snpList[i]->GetrsId())==snpIndex.end()){
            snpIndex[snpList[i]->GetrsId()] = i;
            try{
                snpList[i]->computeVarianceExplainedChi(caseSize, controlSize, prevalence, isPvalue);;
                snpList[i]->m_regionFlag.push_back(false);
                if(regionList->GetnumRegion() != 0){
                    std::vector<bool> padding(regionList->GetnumRegion(), false);
                    snpList[i]->m_regionFlag.insert(snpList[i]->m_regionFlag.end(), padding.begin(), padding.end());
                }
                for(size_t j = 0; j < regionList->GetnumRegion(); ++j){
                    for(unsigned k = regionIncrementationIndex.at(j); k < regionList->GetintervalSize(j); ++k){
                        if(regionList->Getchr(j,k).compare(snpList[i]->Getchr())==0 &&
                           regionList->Getstart(j,k) <= snpList[i]->Getbp() &&
                           regionList->Getend(j,k) >= snpList[i]->Getbp()){
                            regionIncrementationIndex.at(j) = k;
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


void Snp::generateSnpIndex(std::map<std::string, size_t> &snpIndex, std::vector<Snp*> &snpList, std::string genotypeFileName, bool keep){
    //If we do keep the ambiguous SNPs, we need to make sure it is the correct orientation
    //not the most sophisticated method though.
    std::string bimFileName = genotypeFileName;
    bimFileName.append(".bim");
    std::ifstream bimFile;
    bimFile.open(bimFileName.c_str());
    if(!bimFile.is_open()){
        std::string message = "Cannot open bim file: ";
        message.append(bimFileName);
        throw message;
    }
    std::string line;
    while(std::getline(bimFile,line)){
        line =usefulTools::trim(line);
    }
    bimFile.close();
	size_t duplicate = 0;
	for(size_t i = 0; i < snpList.size(); ++i){
        if(snpIndex.find(snpList[i]->GetrsId())==snpIndex.end()){
            snpIndex[snpList[i]->GetrsId()] = i;
        }
        else{
            duplicate++;
        }

    }
    if(duplicate == 0) std::cerr << "There are no duplicated rsID in the p-value file" << std::endl << std::endl;
    else std::cerr <<  "There are a total of " << duplicate << " duplicated rsID(s) in the p-value file" << std::endl << std::endl;
}



void Snp::addDirection(std::map<std::string, size_t> &snpIndex, std::vector<Snp*> &snpList,std::string dirFile){
    if(dirFile.empty()){
        throw "Should not happen as we direction file was not requested";
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
                if(snpIndex.find(token.at(0)) != snpIndex.end()){
                    size_t refId = snpIndex.at(token.at(0));
                    snpList.at(refId)->Setsign(atoi(line.c_str()));
                }
                else{
                    throw "The direction file contains Snps not found in p-value file. Most likely they are not matched. Please check your input!";
                }
            }
        }
    }
    direction.close();
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


void Snp::computeVarianceExplainedChi(bool isPvalue){
    std::cerr << "not supported yet!" << std::endl;
    exit(-1);
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
    m_sampleSize = 20000;
    (*m_beta) = (*m_beta)*(*m_beta);
    (*m_sqrtChiSq)= sqrt((*m_beta))*m_sign;
	(*m_beta) = ((*m_beta)-1)/(m_sampleSize-2.0+(*m_beta));
}



void Snp::setFlag(size_t index, bool value){
    m_regionFlag.at(index) = value;
}
void Snp::addFlag(bool value){
    m_regionFlag.push_back(value);
}

void Snp::cleanSnp(std::vector<Snp*> &snpList){
    for(size_t i= 0; i < snpList.size(); ++i){
        delete snpList.at(i);
    }
    snpList.clear();
}

bool Snp::Concordant(std::string chr, size_t bp, std::string rsId) const{
    return chr.compare(m_chr) ==0 && bp==m_bp && rsId.compare(m_rs) == 0;
}

bool Snp::ambiguousAllele(const std::string refAllele, const std::string altAllele){
    if( (refAllele.compare("A")==0 || refAllele.compare("a")==0) &&
        (altAllele.compare("T")==0 || refAllele.compare("t")==0)) return true;
    if( (refAllele.compare("T")==0 || refAllele.compare("t")==0) &&
        (altAllele.compare("A")==0 || refAllele.compare("a")==0)) return true;
    if( (refAllele.compare("G")==0 || refAllele.compare("G")==0) &&
        (altAllele.compare("C")==0 || refAllele.compare("c")==0)) return true;
    if( (refAllele.compare("C")==0 || refAllele.compare("c")==0) &&
        (altAllele.compare("G")==0 || refAllele.compare("g")==0)) return true;
    return false;
}






