#include "snpestimation.h"

SnpEstimation::SnpEstimation(GenotypeFileHandler *genotypeFileHandler, SnpIndex *snpIndex, std::vector<Snp*> *snpList, size_t thread, double maf, bool correction):m_genotypeFileHandler(genotypeFileHandler), m_snpIndex(snpIndex), m_snpList(snpList), m_thread(thread), m_maf(maf), m_correction(correction){};

void SnpEstimation::Estimate(){
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
		std::cerr << "Getting Snps" << std::endl;
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
			std::cerr << "Compute LD" << std::endl;
			ProcessCode linkageProcess = linkageMatrix->Construct(genotype, prevResidual, blockSize, m_correction);
			if(linkageProcess == fatalError){
                exit(-1);
            }
            //Now we can perform the decomposition on the data
            std::cerr << "Decompose" << std::endl;
            decompositionHandler->Decompose(blockSize, snpLoc, genotype, chromosomeStart, chromosomeEnd);
            if(blockSize > genotype.size()) blockSize= genotype.size();
            Genotype::clean(genotype, blockSize);
            size_t removeCount = snpLoc.size() - blockSize;
            for(size_t i = 0; i < removeCount; ++i)	snpLoc.pop_front();


		}
        if(chromosomeStart && !chromosomeEnd){
            chromosomeStart =false;
        }
        else if(chromosomeEnd){
            chromosomeStart = true;
            chromosomeEnd = false;
        }
	}
	m_effective = linkageMatrix->Geteffective();
	m_bug = linkageMatrix->Getbug();//DEBUG
	delete linkageMatrix;
	delete decompositionHandler;
    Genotype::clean(genotype, 0);
}

SnpEstimation::~SnpEstimation()
{
	//dtor
}

void SnpEstimation::Getresult(std::string outputPrefix){
    size_t regionSize = m_snpList->front()->GetregionSize();
    std::vector<double> regionEstimate(regionSize, 0.0);
    size_t nSnp = 0;
    m_snpIndex->init();
    size_t index;
    double totalSum = 0.0;
    while(m_snpIndex->valid()){
		index = m_snpIndex->value();
        double num = (*m_snpList)[index]->Getheritability();
        totalSum+= num;
        if((*m_snpList)[index]->GetFlag(0)){
			nSnp++;
            for(size_t j = 0; j < regionSize;++j){
                if((*m_snpList)[index]->GetFlag(j)){
                    regionEstimate[j]+= num;
                }
            }
        }
        if(!m_snpIndex->next()){
            break;
        }

    }
    std::cerr << "Effective Number: " << (m_effective) << std::endl; //DEBUG
    std::cerr << "Effective Number3: " <<m_bug << std::endl; //DEBUG
	if(outputPrefix.empty()){
        std::cout << "Category\tPositive\tNegative\tVariance" << std::endl;
        std::cout << "With LD\t" << regionEstimate[0] << "\t" << totalSum-regionEstimate[0] << "\t" << (2.0*(nSnp*nSnp/m_effective)+4.0*m_snpList->front()->GetsampleSize()*regionEstimate[0])/(m_snpList->front()->GetsampleSize()*m_snpList->front()->GetsampleSize()*1.0) << std::endl;
        for(size_t i =1; i < regionEstimate.size(); ++i){
            std::cout << Region::regionNames[i] << "\t" << regionEstimate[i] << "\t" << totalSum-regionEstimate[i] << "\tNA" << std::endl;
        }
    }
    else{
		std::string resOutName = outputPrefix+".res";
        std::string resSumName = outputPrefix+".sum";
        std::ofstream resOut, resSum;
        resOut.open(resOutName.c_str());
        resSum.open(resSumName.c_str());


        if(!resOut.is_open()){
            std::cerr << "Cannot open result file: " << resOutName << " for write" << std::endl;
            std::cerr << "Will only provide brief summary output " << std::endl;
        }
        else{
            resOut << "Chr\tLoc\trsID\tOriginal\tBeta\tEstimate\tWithLD" << std::endl;
            m_snpIndex->init();
            size_t index;
            while(m_snpIndex->valid()){
                index = m_snpIndex->value();
                resOut << (*m_snpList)[index]->Getchr() << "\t" << (*m_snpList)[index]->Getbp() << "\t" << (*m_snpList)[index]->GetrsId() << "\t" << (*m_snpList)[index]->Getoriginal()<< "\t" << (*m_snpList)[index]->Getbeta() << "\t" <<  (*m_snpList)[index]->Getheritability() << "\t" << (*m_snpList)[index]->GetFlag(0) << std::endl;
                if(!m_snpIndex->next()){
                    break;
                }
            }
        }
        if(!resSum.is_open()){
            std::cerr << "Cannot open summary file: " << resSumName << " for write" << std::endl;
            std::cerr << "Will display on screen" << std::endl;
            std::cout << "Category\tPositive\tNegative\tVariance" << std::endl;
            std::cout << "With LD\t" << regionEstimate[0] << "\t" << totalSum-regionEstimate[0] << "\t" << (2.0*(nSnp*nSnp/m_effective)+4.0*m_snpList->front()->GetsampleSize()*regionEstimate[0])/(m_snpList->front()->GetsampleSize()*m_snpList->front()->GetsampleSize()*1.0) << std::endl;
			for(size_t i =1; i < regionEstimate.size(); ++i){
				std::cout << Region::regionNames[i] << "\t" << regionEstimate[i] << "\t" << totalSum-regionEstimate[i] << "\tNA" << std::endl;
			}
        }
        else{
            resSum << "Category\tPositive\tNegative\tVariance" << std::endl;
            resSum << "With LD\t" << regionEstimate[0] << "\t" << totalSum-regionEstimate[0] << "\t" << (2.0*(nSnp*nSnp/m_effective)+4.0*m_snpList->front()->GetsampleSize()*regionEstimate[0])/(m_snpList->front()->GetsampleSize()*m_snpList->front()->GetsampleSize()*1.0) << std::endl;
            std::cout << "Category\tPositive\tNegative\tVariance" << std::endl;
            std::cout << "With LD\t" << regionEstimate[0] << "\t" << totalSum-regionEstimate[0] << "\t" << (2.0*(nSnp*nSnp/m_effective)+4.0*m_snpList->front()->GetsampleSize()*regionEstimate[0])/(m_snpList->front()->GetsampleSize()*m_snpList->front()->GetsampleSize()*1.0) << std::endl;
            for(size_t i =1; i < regionEstimate.size(); ++i){
				resSum << Region::regionNames[i] << "\t" << regionEstimate[i] << "\t" << totalSum-regionEstimate[i] << "\tNA" << std::endl;
				std::cout << Region::regionNames[i] << "\t" << regionEstimate[i] << "\t" << totalSum-regionEstimate[i] << "\tNA" << std::endl;
			}
            resSum.close();
        }
    }
}
