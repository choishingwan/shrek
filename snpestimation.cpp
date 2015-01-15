#include "snpestimation.h"

SnpEstimation::SnpEstimation(GenotypeFileHandler *genotypeFileHandler, SnpIndex *snpIndex, std::vector<Snp*> *snpList, size_t blockSize, size_t distance, size_t thread, double maf):m_genotypeFileHandler(genotypeFileHandler), m_snpIndex(snpIndex), m_snpList(snpList), m_blockSize(blockSize), m_distance(distance), m_thread(thread), m_maf(maf){};

void SnpEstimation::performEstimation(){
/**
 *	while(c=function(), c==test)
 *	will first perform function and assign output to c, then return the value of the comparison of c and test (exactly what we want
 */
 /**
  * Things to pay attention to:
  * 1. Snp density are not completely even
  * 2. Need to be-careful with possible duplicated Snps within the LD file
  * 3. What is a "block"? <- most difficult question
  * 4. Need to consider when there is no snps with LD information
  */
	Genotype::SetSampleNum(m_genotypeFileHandler->GetSampleSize());
	ProcessCode process = startProcess;
    std::deque<Genotype*> genotype;
    std::deque<size_t> snpLoc;
	while(process != completed){
		m_genotypeFileHandler->getSnps(genotype, m_distance, 0, "", m_snpList);
        //p-impute
        //remove snps
	}
}

SnpEstimation::~SnpEstimation()
{
	//dtor
}
