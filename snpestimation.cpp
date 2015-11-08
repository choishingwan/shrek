#include "snpestimation.h"

SnpEstimation::SnpEstimation(){};
SnpEstimation::~SnpEstimation(){};

void SnpEstimation::Estimate(GenotypeFileHandler &genotypeFileHandler,const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, const Region &regionInfo, const Command &commander,boost::ptr_vector<Interval> &blockInfo){
    //Declaration

	Genotype::SetsampleNum(genotypeFileHandler.getSampleSize()); //Cannot forget this, otherwise the programme will crash due to not knowing the sample size
    //TODO: Make this more trivial
    size_t genotypeResidual = 0, startBlockIndex = 0;; //This is index of the block to read SNPs
    size_t previousLeftOvers=0;
    boost::ptr_deque<Genotype> genotype;
    std::deque<size_t> snpLoc; //Store the SNP index (for snpList)
    std::deque<size_t> ldLoc; //Store the Genotype Index, for use with blockInfo which gives bound of this index
    std::deque<size_t> blockLoc;
    bool chromosomeStart = true;
    bool chromosomeEnd = false;
    bool correction = commander.ldCorrect();
    Linkage linkageMatrix(commander.getThread());
    //bool direction = commander.hasDir();
    bool direction = false;
    Decomposition decomposition(commander.getThread(), direction);

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
        //All getSnps need is an end of reading
        //std::cerr << "Start getting SNPs: " << previousLeftOvers << "\t" << genotypeResidual << std::endl;
        genotypeFileHandler.getSnps(genotype, snpLoc, ldLoc, chromosomeStart, chromosomeEnd, genotypeResidual, blockInfo);
        //std::cerr << "Now got: " << genotype.size() << std::endl;
        /** The interval of blockInfo should be [) not [] **/
        /** The genotypeResidual is ALWAYS pointing to the NEXT BLOCK **/
        //build linkage
        linkageMatrix.Initialize(genotype, previousLeftOvers);
        linkageMatrix.Construct(genotype, startBlockIndex, previousLeftOvers, blockInfo, correction, ldLoc);
        //linkageMatrix.print();
        //exit(-1);
        decomposition.run(linkageMatrix, startBlockIndex, previousLeftOvers, blockInfo, ldLoc, snpLoc, snpList, chromosomeStart);
        //std::cerr << "Decomposed" << std::endl;
        //Now we have the LD matrix we will want to decompose it
        //If it is not chromosome end, we remove everything except the last two blocks
        if(chromosomeEnd){
            chromosomeEnd = false;
            chromosomeStart = true;
            //When end, it means all things in this chromosome is over, and we should not use the current block

            genotype.clear();
            ldLoc.clear();
            snpLoc.clear();
        }
        else{
            chromosomeStart = false;
            //Find everything except last two
            //genotypeResidual is the last block+1
            size_t startKeep = blockInfo[genotypeResidual-2].getStart();
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

}

void SnpEstimation::getResult(const Command &commander, const Region &region, const std::map<std::string, size_t> &snpIndex, const boost::ptr_vector<Snp> &snpList){
    size_t regionSize = region.getNumRegion();
    std::vector<double> regionEstimate(regionSize, 0.0);
    std::vector<double> regionEffect(regionSize, 0.0);
    std::map<std::string, size_t>::const_iterator iter;
    double totalSum = 0.0;
    size_t index=0;
    size_t sampleSize = 1;
    for(iter= snpIndex.begin(); iter != snpIndex.end(); ++iter){
		index =iter->second;
        double num = snpList.at(index).getHeritability();
        double eff =snpList.at(index).getEffectiveNumber();
        totalSum+= num;
        if(snpList.at(index).getFlag(0)){
            for(size_t j = 0; j < regionSize;++j){
                if(snpList.at(index).getFlag(j)){
                    regionEstimate[j]+= num;
                    regionEffect[j]+= eff;
                }
            }
        }
    }
    double adjustment = 1.0;
    if(commander.quantitative()){
        adjustment = commander.getExtreme();
        if(commander.sampleSizeProvided()) sampleSize =  commander.getSampleSize();
        else sampleSize = Snp::getMaxSample();
    }
    else if(commander.caseControl()){
        sampleSize = commander.getCaseSize()+commander.getControlSize();
        double portionCase = (double)(commander.getCaseSize()) / (double)(sampleSize);
        double i2 = usefulTools::dnorm(usefulTools::qnorm(commander.getPrevalence()))/(commander.getPrevalence());
        i2 = i2*i2;
        adjustment = ((1.0-commander.getPrevalence())*(1.0-commander.getPrevalence()))/(i2*portionCase*(1-portionCase));
    }
    totalSum *= adjustment;
    for(size_t i = 0; i <regionEstimate.size(); ++i){
        //can't remember why I do the adjustment^2 and adjustment separately, something to do with our old implementation I guess
        regionEffect[i] = 2.0*(regionEffect[i]*adjustment*adjustment+2.0*adjustment*regionEstimate[i]*sampleSize)/((double)sampleSize*sampleSize);
        regionEstimate[i] *= adjustment;
    }

    std::cout << "Category\tPositive\tNegative\tVariance" << std::endl;
    for(size_t i =0; i < regionEstimate.size(); ++i){
        std::cout << region.getName(i) << "\t" << regionEstimate[i] << "\t" << totalSum-regionEstimate[i] << "\t" << regionEffect[i]  << std::endl;
    }
	if(!commander.getOutputPrefix().empty()){
        std::string resOutName = commander.getOutputPrefix()+".res";
        std::string resSumName = commander.getOutputPrefix()+".sum";
        std::ofstream resOut, resSum;
        resOut.open(resOutName.c_str());
        resSum.open(resSumName.c_str());
        if(!resOut.is_open()){
            std::cerr << "Cannot open result file: " << resOutName << " for write" << std::endl;
            std::cerr << "Will only provide brief summary output " << std::endl;
        }
        else{
			//We want the output in sorted format. So we will have to do extra works just to make sure the result is sorted
            resOut << "Chr\tLoc\trsID\tOriginal\tBeta\tEstimate\tEffectiveNumber\tldsc" << std::endl;
            size_t index;
            std::vector<size_t> resultSnpIndex;
            for(iter = snpIndex.begin(); iter != snpIndex.end(); ++iter){
                resultSnpIndex.push_back(index = iter->second);
            }
            std::sort(resultSnpIndex.begin(), resultSnpIndex.end());
            for(size_t i =0; i < resultSnpIndex.size(); ++i){
                size_t index =resultSnpIndex[i];
				resOut << snpList[index].getChr() << "\t" << snpList[index].getBp() << "\t" << snpList[index].getRs() << "\t" << snpList[index].getOriginal()<< "\t" << snpList[index].getBeta() << "\t" <<  snpList[index].getHeritability() << "\t" <<  snpList[index].getEffectiveNumber() <<"\t" << snpList[index].getLDScore()<< std::endl;
            }
        }

        if(!resSum.is_open()){
            std::cerr << "Cannot open summary file: " << resSumName << " for write" << std::endl;
        }
        else{
            resSum << "Category\tPositive\tNegative\tVariance" << std::endl;
			for(size_t i =0; i < regionEstimate.size(); ++i){
				resSum << region.getName(i) << "\t" << regionEstimate[i] << "\t" << totalSum-regionEstimate[i] << "\t" << regionEffect[i]  << std::endl;
			}
            resSum.close();
        }
    }

}
