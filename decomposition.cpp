#include "decomposition.h"

Decomposition::Decomposition(SnpIndex *snpIndex, std::vector<Snp*> *snpList, Linkage *linkageMatrix, size_t thread):m_snpIndex(snpIndex), m_snpList(snpList), m_linkage(linkageMatrix), m_thread(thread){}

Decomposition::~Decomposition()
{
	//dtor
}

ProcessCode Decomposition::Decompose(const size_t &blockSize, std::deque<size_t> &snpLoc, std::deque<Genotype*> &genotype, bool chromosomeStart, bool chromosomeEnd){
    //First, build the vector for beta
    if(genotype.size() == 0){
		return continueProcess;
    }

    size_t processSize = snpLoc.size();
    if(!chromosomeEnd){ //We have not reached the end, so we need to remove the last two blocks from the processing
        processSize = blockSize/3*m_thread+2*(blockSize/3);
    }
    //std::cerr << "Chromosome is end: " << chromosomeEnd << std::endl;
    //std::cerr << "snpLoc size: " << snpLoc.size() << std::endl;
    Eigen::VectorXd betaEstimate = Eigen::VectorXd::Zero(processSize);
    Eigen::VectorXd signValue = Eigen::VectorXd::Constant(processSize, 1.0);
	for(size_t i=0;i < processSize; ++i){
        if(snpLoc[i] < 0){
			std::cerr << "ERROR! Undefined behaviour. The location index is negative!" << std::endl;
			return fatalError;
        }
        else{
			betaEstimate(i) = (*m_snpList)[snpLoc[i]]->Getbeta();
			signValue(i) = (*m_snpList)[snpLoc[i]]->Getsign();
        }
    }
    //Now we can perform decomposition using the beta and Linkage

	size_t currentBlockSize = blockSize;
	if(currentBlockSize > processSize){
        currentBlockSize = processSize;
	}
    std::vector<size_t> startLoc;
	size_t stepSize = currentBlockSize/3;
    if(stepSize == 0 || currentBlockSize==processSize){ //no multithreading is required
		//The Block size is only 3. so We can finish it anyway
		Eigen::MatrixXd variance = Eigen::MatrixXd::Zero(processSize, processSize);
		Eigen::VectorXd result = m_linkage->solveChi(0, processSize, &betaEstimate, &signValue, &variance, Snp::m_maxSampleSize);
		size_t copyStart = 0;
		if(!chromosomeStart) copyStart = blockSize/3;
		std::cerr << "Region size: " << Region::regionVariance.size() << std::endl;
        for(size_t i=copyStart; i < processSize; ++i){
			(*m_snpList)[snpLoc[i]]->Setheritability(result(i));
			(*m_snpList)[snpLoc[i]]->Setvariance(variance(i,i)); //The diagonal of the matrix contains the per snp variance
			for(size_t j = 0; j < processSize; ++j){ //For this snp, go through all the snp partners
				double covariance = variance(i,j);
				for(size_t regionIndex = 0; regionIndex < Region::regionVariance.size(); ++regionIndex){
					if((*m_snpList)[snpLoc[i]]->GetFlag(regionIndex)&&(*m_snpList)[snpLoc[j]]->GetFlag(regionIndex)){
						Region::regionVariance[regionIndex]+= covariance;
					}
				}
			}
        }

	}
	else if(stepSize > 0){
			//std::cerr << "Multithread " << std::endl;
		//Normal situation
		for(size_t i = 0; i < processSize; i+= stepSize){
			if(i+currentBlockSize <= processSize){
				startLoc.push_back(i);
			}
		}
		if(!startLoc.empty() && processSize - startLoc.back() < stepSize){
			startLoc.pop_back();
		}
		std::vector<pthread_t*> threadList;
		std::vector<DecompositionThread* > garbageCollection;
		if(m_thread < startLoc.size() && !chromosomeEnd){ //we should be only using the first part, so there shouldn't be any problem
			std::cerr << "Undefined behaviour! Don't worry, this should be the problem of the programmer instead of the user..." << std::endl;
			std::cerr << startLoc.size() << "\t" << m_thread << std::endl;
			return fatalError;
		}
		else if(m_thread < startLoc.size()){
            //std::cerr << "Extra" << std::endl;
            //We will first do the appropriate amount, then do the remaining;
            for(size_t i =0; i < m_thread; ++i){
				garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate, &signValue, m_linkage, &snpLoc, m_snpList, chromosomeStart, false));
				if(i == startLoc.size()-1) garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate,&signValue, m_linkage, &snpLoc, m_snpList, chromosomeStart, true));
				else garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate,&signValue, m_linkage, &snpLoc, m_snpList, chromosomeStart, false));
				pthread_t *thread1 = new pthread_t();
				threadList.push_back(thread1);
				int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection.back());
				if(threadStatus != 0){
					std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
					return fatalError;
				}
            }
			for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
			for(size_t i = 0; i < threadList.size(); ++i) delete threadList[i];
            for(size_t i = 0; i < garbageCollection.size(); ++i) delete garbageCollection[i];
			garbageCollection.clear();
            threadList.clear();
            //Now we do the extra job
			for(size_t i=m_thread; i < startLoc.size(); ++i){
				if(i == startLoc.size()-1) garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate,&signValue, m_linkage, &snpLoc, m_snpList, chromosomeStart, true));
				else garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate,&signValue, m_linkage, &snpLoc, m_snpList, chromosomeStart, false));
				pthread_t *thread1 = new pthread_t();
				threadList.push_back(thread1);
				int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection.back());
				if(threadStatus != 0){
					std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
					return fatalError;
				}
            }
            for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
			for(size_t i = 0; i < threadList.size(); ++i) delete threadList[i];
            for(size_t i = 0; i < garbageCollection.size(); ++i) delete garbageCollection[i];
            threadList.clear();
            garbageCollection.clear();

		}
		else{
			for(size_t i = 0; i < startLoc.size(); ++i){
				if(i == startLoc.size()-1) garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate,&signValue,m_linkage, &snpLoc, m_snpList, chromosomeStart, true));
				else garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate,&signValue, m_linkage, &snpLoc, m_snpList, chromosomeStart, false));
				pthread_t *thread1 = new pthread_t();
				threadList.push_back(thread1);
				int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection.back());
				if(threadStatus != 0){
					std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
					return fatalError;
				}
			}
			for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
			for(size_t i = 0; i < threadList.size(); ++i) delete threadList[i];
            for(size_t i = 0; i < garbageCollection.size(); ++i) delete garbageCollection[i];
            threadList.clear();
            garbageCollection.clear();
		}
	}
	else{
        std::cerr << "Undefined behaviour! Step size should never be negative as block size is positive" << std::endl;
        return fatalError;
	}

	return completed;
}

