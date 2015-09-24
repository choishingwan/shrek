#include "snpestimation.h"

SnpEstimation::SnpEstimation(){};
SnpEstimation::~SnpEstimation(){};

void SnpEstimation::Estimate(GenotypeFileHandler *genotypeFileHandler,const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, Region* regionInfo, Command *commander,boost::ptr_vector<Interval> &blockInfo){
    //Declaration

	Genotype::SetsampleNum(genotypeFileHandler->getSampleSize()); //Cannot forget this, otherwise the programme will crash due to not knowing the sample size
    //TODO: Make this more trivial
    size_t genotypeResidual = 0; //This is index of the block to read SNPs
    size_t remainLd = 0; //This is the index of the block to construct LD
    boost::ptr_deque<Genotype> genotype;
    std::deque<size_t> snpLoc;
    std::deque<size_t> ldLoc;
    bool chromosomeStart = true;
    bool chromosomeEnd = false;
    bool correction = commander->ldCorrect();
    Linkage *linkageMatrix = nullptr;
    linkageMatrix = new Linkage(commander->getThread());
    //Start processing (will need to wrap it with a while loop
    /*
    1. Get SNPs
    2. Build Linkage
    3. Decomposition
    */
    //The following
    while(genotypeResidual != blockInfo.size()){ //When we reaches blockInfo.size, it means we have finish all work
        genotypeFileHandler->getSnps(genotype, snpLoc, ldLoc, snpList, chromosomeStart, chromosomeEnd, genotypeResidual, blockInfo);


        //build linkage
        //linkageMatrix->Initialize(genotype, remainLd, blockInfo);
        //linkageMatrix->Construct(genotype, genotypeResidual, remainLd, blockInfo, correction, ldLoc);
        //Now we have the LD matrix we will want to decompose it
    }


}
