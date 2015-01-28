#include "snpestimation.h"

SnpEstimation::SnpEstimation(GenotypeFileHandler *genotypeFileHandler, SnpIndex *snpIndex, std::vector<Snp*> *snpList, size_t thread, double maf, bool correction):m_genotypeFileHandler(genotypeFileHandler), m_snpIndex(snpIndex), m_snpList(snpList), m_thread(thread), m_maf(maf), m_correction(correction){};

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
	Genotype::SetsampleNum(m_genotypeFileHandler->GetsampleSize());
	ProcessCode process = startProcess;
    std::deque<Genotype*> genotype;
    std::deque<size_t> snpLoc;
    bool chromosomeStart= true;
    bool chromosomeEnd = false;
    size_t prevResidual;
    size_t blockSize;
    Linkage *linkageMatrix = new Linkage(m_thread);
    Decomposition *decompositionHandler = new Decomposition( m_snpIndex, m_snpList, linkageMatrix, m_thread);
	while(process != completed && process != fatalError){
		process = m_genotypeFileHandler->getSnps(genotype, snpLoc, *m_snpList, chromosomeStart, chromosomeEnd, m_maf,prevResidual, blockSize);
		if(process == fatalError){
            exit(-1);
		}
		if(process == completed && prevResidual==genotype.size()){
			//Nothing was updated
			std::cerr << "completed" << std::endl;
		}
		else{
			//Now calculate the LD matrix
			linkageMatrix->Construct(genotype, prevResidual, blockSize, m_correction);
            //Now we can perform the decomposition on the data

		}
        //p-impute
        //remove snps
	}
	delete linkageMatrix;
	delete decompositionHandler;
}

SnpEstimation::~SnpEstimation()
{
	//dtor
}
