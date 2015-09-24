#include "snpestimation.h"

SnpEstimation::SnpEstimation(){};

void SnpEstimation::Estimate(GenotypeFileHandler *genotypeFileHandler,const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, Region* regionInfo, Command *commander,boost::ptr_vector<Interval> &blockInfo){
    //Declaration
    size_t prevResidual = 0; //This is the index of the last finished block
    boost::ptr_deque<Genotype> genotype;
    std::deque<size_t> snpLoc;
    bool chromosomeStart = true;
    bool chromosomeEnd = false;
    //Start processing (will need to wrap it with a while loop
    /*
    1. Get SNPs
    2. Build Linkage
    3. Decomposition
    */
    //The following
    genotypeFileHandler->getSnps(genotype, snpLoc, snpList, chromosomeStart, chromosomeEnd, prevResidual, blockInfo);
    //build linkage

}

	Genotype::SetsampleNum(m_genotypeFileHandler->getSampleSize());
    boost::ptr_deque<Genotype> genotype;
    std::deque<size_t> snpLoc;

    bool chromosomeStart= true;
    bool chromosomeEnd = false;

    Linkage *linkageMatrix = new Linkage();
    linkageMatrix->setSnpList(m_snpList);
    linkageMatrix->setSnpLoc(&snpLoc);
    linkageMatrix->setThread(m_thread);

    Decomposition *decompositionHandler = new Decomposition( m_snpList, linkageMatrix, m_thread,m_regionInfo);
	while(process != completed && process != fatalError){
		/** Will have terrible problem if the input is corrupted */
        SnpEstimation::loadbar(numProcessed,totalNum);
        //Work one by one first.
		process = m_genotypeFileHandler->getSnps(genotype, snpLoc, m_snpList, chromosomeStart, chromosomeEnd, m_maf,prevResidual, blockSize);

		if(process == completed && !chromosomeEnd){
			m_regionInfo->Debuffer();
		}
		else{
			//Now calculate the LD matrix
			linkageMatrix->Initialize(genotype, prevResidual, blockSize);
			linkageMatrix->Construct(genotype, prevResidual, blockSize, m_correction);

            m_regionInfo->CleanBuffer();
            ProcessCode decomposeProcess = decompositionHandler->Decompose(blockSize, snpLoc, genotype, chromosomeStart, chromosomeEnd);

			if(decomposeProcess == fatalError) throw "Fatal error with Decomposition";
            numProcessed+= genotype.size(); //Finished the LD construction
            if(!chromosomeEnd){
				if(blockSize > genotype.size()) throw "When block size is bigger than the number of genotype, it must be the end of chromosome";
				size_t retain = blockSize/3*2;
				Genotype::clean(genotype, retain);
				size_t removeCount = snpLoc.size() - retain;
				for(size_t i = 0; i < removeCount; ++i)	snpLoc.pop_front();
				numProcessed-=retain;
            }
            else{
                m_regionInfo->Debuffer();
                Genotype::clean(genotype,0);
                snpLoc.clear();
            }

            SnpEstimation::loadbar(numProcessed,totalNum);
		}
        if(chromosomeStart && !chromosomeEnd){
            chromosomeStart =false;
        }
        else if(chromosomeEnd){
            m_regionInfo->Debuffer();
            chromosomeStart = true;
            chromosomeEnd = false;
        }
	}
/*
	std::cerr << "Here" << std::endl;
	std::cout << DecompositionThread::checking << std::endl;
    linkageMatrix->print();
	exit(-1);
*/
    //std::cout << Linkage::m_testing << std::endl;
	m_regionInfo->Debuffer();
    SnpEstimation::loadbar(totalNum, totalNum);
	std::cerr << std::endl;
	delete linkageMatrix;
	delete decompositionHandler;
    Genotype::clean(genotype, 0);
}
