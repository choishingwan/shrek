#include "snp.h"

Snp::Snp(std::string chr, std::string rs, size_t bp, size_t sampleSize, double original, double beta):m_chr(chr), m_rs(rs), m_bp(bp), m_sampleSize(sampleSize), m_original(original), m_beta(beta){
	m_heritability = 0.0;
}

std::string Snp::Getchr() const { return m_chr; }
std::string Snp::GetrsId() const { return m_rs; }
size_t Snp::Getbp() const { return m_bp; }
size_t Snp::GetsampleSize() const { return m_sampleSize; }
size_t Snp::GetregionSize() const {return m_regionFlag.size(); }
double Snp::Getoriginal() const { return m_original; }
double Snp::Getbeta() const { return m_beta; }
double Snp::Getheritability() const { return m_heritability; }
void Snp::Setheritability(double heritability ) { m_heritability = heritability; }
bool Snp::GetFlag(size_t index) const {
	if(index>= m_regionFlag.size()) return false;
	return m_regionFlag[index];
}

Snp::~Snp()
{
	//dtor
}


void Snp::generateSnpList(std::vector<Snp*> &snpList, const std::string &pvalueFile, const size_t index, const size_t sampleSize, const size_t rsIndex, const size_t bpIndex, const size_t chrIndex, const size_t sampleIndex, bool sampleSizeProvided){
	std::ifstream pValue;
    pValue.open(pvalueFile.c_str());
    if(!pValue.is_open()){
		std::cerr << "Cannot open the p-value file: " << pvalueFile << std::endl;
		exit(-1);
    }
    std::string line;
    std::getline(pValue, line); //First line is header
    while(std::getline(pValue, line)){
		line = usefulTools::trim(line);
        if(!line.empty()){
			std::vector<std::string> token;
			usefulTools::tokenizer(line, " \t", &token);
            if(token.size() > index && token.size() > rsIndex && token.size() > bpIndex && token.size() > chrIndex){
				if(sampleSizeProvided && token.size() <= sampleIndex){
					//Skip, as the sample index will be out of bound
				}
				else{
					std::string chr = token[chrIndex];
                    size_t bp = atoi(token[bpIndex].c_str());
                    std::string rsId = token[rsIndex];
                    size_t sizeOfSample = sampleSize;
                    if(!sampleSizeProvided) sizeOfSample = atoi(token[sampleIndex].c_str());
                    double predictedBeta = atof(token[index].c_str());
					if(!std::isfinite(predictedBeta)){
                        predictedBeta=0;
                        std::cerr << rsId << " does not have finite input, will set it to zero" << std::endl;
                    }
                    snpList.push_back(new Snp(chr, rsId, bp, sizeOfSample, predictedBeta, predictedBeta));
				}
            }
        }
    }
    pValue.close();
    std::sort(snpList.begin(), snpList.end(), Snp::sortSnp);
    snpList.erase( unique( snpList.begin(), snpList.end() ), snpList.end() );
    std::cerr << "There are a total of " << snpList.size() << " Snps in the input" << std::endl;
}

bool Snp::sortSnp (Snp* i, Snp* j){
    if(i->Getchr().compare(j->Getchr()) == 0)
		if(i->Getbp() == j->Getbp())
			return i->GetrsId().compare(j->GetrsId()) < 0;
		else
			return i->Getbp() < j->Getbp();
	else return (i->Getchr().compare(j->Getchr()) < 0);
}

void Snp::generateSnpIndex(SnpIndex *snpIndex, std::vector<Snp*> &snpList, std::vector<std::vector<Region*> > &regionList, bool isPvalue){
	std::vector<size_t> regionIncrementationIndex(regionList.size(), 0);
	size_t duplicate = 0;
	for(size_t i = 0; i < snpList.size(); ++i){
        if(!snpIndex->find(snpList[i]->GetrsId())){
            snpIndex->set(snpList[i]->GetrsId(), i);
            snpList[i]->computeVarianceExplained(isPvalue);
            snpList[i]->m_regionFlag.push_back(false);

            if(regionList.size() != 0){
                std::vector<bool> padding(regionList.size(), false);
                snpList[i]->m_regionFlag.insert(snpList[i]->m_regionFlag.end(), padding.begin(), padding.end());
            }
            for(size_t j = 0; j < regionList.size(); ++j){
                for(unsigned k = regionIncrementationIndex[j]; k < regionList[j].size(); ++k){
                    if(regionList[j][k]->Getchr().compare(snpList[i]->Getchr())==0 &&
                       regionList[j][k]->Getstart() <= snpList[i]->Getbp() &&
                       regionList[j][k]->Getend() >= snpList[i]->Getbp()){
                        //std::cerr << snpList[i]->GetrsId() << "\tWithin region" << std::endl;
                        regionIncrementationIndex[j] = k;
                        snpList[i]->m_regionFlag[j+1]=true;
                        break;
                       }
                }
            }
        }
        else{
            duplicate++;
        }

    }
    if(duplicate == 0) std::cerr << "There are no duplicated rsID in the p-value file" << std::endl << std::endl;
    else std::cerr <<  "There are a total of " << duplicate << " duplicated rsID(s) in the p-value file" << std::endl << std::endl;
}

void Snp::generateSnpIndex(SnpIndex *snpIndex, std::vector<Snp*> &snpList, const size_t &caseSize, const size_t &controlSize, const double &prevalence, std::vector<std::vector<Region*> > &regionList, bool isPvalue){
	std::vector<size_t> regionIncrementationIndex(regionList.size(), 0);
	size_t duplicate = 0;
	for(size_t i = 0; i < snpList.size(); ++i){
        if(!snpIndex->find(snpList[i]->GetrsId())){
            snpIndex->set(snpList[i]->GetrsId(), i);
            snpList[i]->computeVarianceExplained(caseSize, controlSize, prevalence, isPvalue);;
            snpList[i]->m_regionFlag.push_back(false);

            if(regionList.size() != 0){
                std::vector<bool> padding(regionList.size(), false);
                snpList[i]->m_regionFlag.insert(snpList[i]->m_regionFlag.end(), padding.begin(), padding.end());
            }
            for(size_t j = 0; j < regionList.size(); ++j){
                for(unsigned k = regionIncrementationIndex[j]; k < regionList[j].size(); ++k){
                    if(regionList[j][k]->Getchr().compare(snpList[i]->Getchr())==0 &&
                       regionList[j][k]->Getstart() <= snpList[i]->Getbp() &&
                       regionList[j][k]->Getend() >= snpList[i]->Getbp()){
                        //std::cerr << snpList[i]->GetrsId() << "\tWithin region" << std::endl;
                        regionIncrementationIndex[j] = k;
                        snpList[i]->m_regionFlag[j+1]=true;
                        break;
                       }
                }
            }
        }
        else{
            duplicate++;
        }

    }
    std::cerr <<  "There are a total of " << duplicate << " duplicated rsID(s)" << std::endl << std::endl;

}


void Snp::computeVarianceExplained(bool isPvalue){
    if(isPvalue){
        m_beta = usefulTools::qnorm(1.0-((m_original+0.0)/2.0));
        if(!std::isfinite(m_beta)) m_beta = usefulTools::qnorm(((m_original+0.0)/2.0));
    }
    m_beta = m_beta*m_beta;
	m_beta = (m_beta/(m_sampleSize-2.0+m_beta))-1.0/(m_sampleSize-1.0);

}

void Snp::computeVarianceExplained(const size_t &caseSize, const size_t &controlSize, const double &prevalence, bool isPvalue){
    if(isPvalue){
        m_beta = usefulTools::qnorm(1.0-((m_original+0.0)/2.0));
        if(!std::isfinite(m_beta)) m_beta =usefulTools::qnorm(((m_original+0.0)/2.0));
        m_beta = m_beta*m_beta;
    }
    double ncp = (m_beta -1.0);
	int totalSampleSize = caseSize + controlSize;
	double portionCase = (caseSize+0.0) / (totalSampleSize+0.0);
	double i2 = usefulTools::dnorm(usefulTools::qnorm(prevalence))/(prevalence);
	i2 = i2*i2;
	m_sampleSize =((1-prevalence)*(1-prevalence))/(i2*portionCase*(1-portionCase)*(totalSampleSize));
	m_beta = m_sampleSize*ncp;

}

void Snp::setFlag(size_t index, bool value){
	if(index >= m_regionFlag.size()){
		std::cerr << "Error: region flag out of bound!" << std::endl;
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