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
    Eigen::VectorXd betaEstimate = Eigen::VectorXd::Zero(m_linkage->rows());
    for(size_t i=0;i < snpLoc.size(); ++i){
        if(snpLoc[i] < 0){
			std::cerr << "ERROR! Undefined behaviour. The location index is negative!" << std::endl;
			return fatalError;
        }
        else betaEstimate(i) = (*m_snpList)[snpLoc[i]]->Getbeta();
    }

    //Now we can perform decomposition using the beta and Linkage

	size_t currentBlockSize = blockSize;
	if(currentBlockSize > genotype.size()){
        currentBlockSize = genotype.size();
	}
    std::vector<size_t> startLoc;
	size_t stepSize = currentBlockSize/3;
	if(stepSize == 0 || currentBlockSize==genotype.size()){ //no multithreading is required
		//The Block size is only 3. so We can finish it anyway
		Eigen::VectorXd result = m_linkage->solve(0, currentBlockSize, &betaEstimate);
        for(size_t i; i < snpLoc.size(); ++i){
			(*m_snpList)[snpLoc[i]]->Setheritability(result(i));
        }

	}
	else if(stepSize > 0){
		//Normal situation
		for(size_t i = 0; i < snpLoc.size(); i+= stepSize){
			if(i+currentBlockSize <= snpLoc.size()){
				startLoc.push_back(i);
			}
		}
		if(!startLoc.empty() && snpLoc.size() - startLoc[startLoc.size()-1] < stepSize){
			startLoc.pop_back();
		}
		std::vector<pthread_t*> threadList;
		std::vector<DecompositionThread* > garbageCollection;
		if(m_thread < startLoc.size()){
			std::cerr << "Undefined behaviour! Don't worry, this should be the problem of the programmer instead of the user..." << std::endl;
			std::cerr << startLoc.size() << "\t" << m_thread << std::endl;
			return fatalError;
		}
		for(size_t i = 0; i < startLoc.size(); ++i){
			if(i == startLoc.size()-1){
				garbageCollection.push_back(new DecompositionThread(startLoc[i], genotype.size()-startLoc[i], &betaEstimate, m_linkage, &snpLoc, m_snpList, chromosomeStart, chromosomeEnd));
			}
			else{
				garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate, m_linkage, &snpLoc, m_snpList, chromosomeStart, chromosomeEnd));
			}
			pthread_t *thread1 = new pthread_t();
			threadList.push_back(thread1);
			int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection.back());
			if(threadStatus != 0){
				std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
				return fatalError;
			}
		}
		for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter){
			pthread_join(*threadList[threadIter], NULL);
		}
		for(size_t i = 0; i < startLoc.size(); ++i){
			delete garbageCollection[i];
			delete threadList[i];
		}
	}
	else{
        std::cerr << "Undefined behaviour! Step size should never be negative as block size is positive" << std::endl;
        return fatalError;
	}

	return completed;
}

