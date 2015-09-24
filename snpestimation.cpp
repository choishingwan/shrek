#include "snpestimation.h"

SnpEstimation::SnpEstimation(){};
SnpEstimation::~SnpEstimation(){};

void SnpEstimation::Estimate(GenotypeFileHandler *genotypeFileHandler,const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, Region* regionInfo, Command *commander,boost::ptr_vector<Interval> &blockInfo){
    //Declaration

	Genotype::SetsampleNum(genotypeFileHandler->getSampleSize()); //Cannot forget this, otherwise the programme will crash due to not knowing the sample size
    //TODO: Make this more trivial
    size_t genotypeResidual = 0, startBlockIndex = 0;; //This is index of the block to read SNPs
    size_t remainLd = 0; //This is the index of the block to construct LD
    size_t previousLeftOvers=0;
    boost::ptr_deque<Genotype> genotype;
    std::deque<size_t> snpLoc; //Store the SNP index (for snpList)
    std::deque<size_t> ldLoc; //Store the Genotype Index, for use with blockInfo which gives bound of this index
    std::deque<size_t> blockLoc;
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
        previousLeftOvers=genotype.size();
        startBlockIndex = genotypeResidual;
        genotypeFileHandler->getSnps(genotype, snpLoc, ldLoc, snpList, chromosomeStart, chromosomeEnd, genotypeResidual, blockInfo);
        /** The interval of blockInfo should be [] not [) **/
        //build linkage
        linkageMatrix->Initialize(genotype, previousLeftOvers);
        linkageMatrix->Construct(genotype, startBlockIndex, previousLeftOvers, blockInfo, correction, ldLoc);
        //Now we have the LD matrix we will want to decompose it

        //If it is not chromosome end, we remove everything except the last two blocks
        if(chromosomeEnd){
            chromosomeEnd = false;
            chromosomeStart = true;
            genotype.clear();
            ldLoc.clear();
            snpLoc.clear();
        }
        else{
            chromosomeStart = false;
            //Find everything except last two
            //genotypeResidual is the last block+1
            size_t startKeep = blockInfo[genotypeResidual-2].getStart();
            size_t endKeep = blockInfo[genotypeResidual-1].getEnd();
            size_t i = ldLoc.size();
            for(;i >0; --i){
                if(ldLoc[i-1] < startKeep){
                    break;
                }
            }
            ldLoc.erase (ldLoc.begin(),ldLoc.begin()+i);
            snpLoc.erase (snpLoc.begin(),snpLoc.begin()+i);
            genotype.erase (genotype.begin(),genotype.begin()+i);
        }

    }
    linkageMatrix->print();
}


