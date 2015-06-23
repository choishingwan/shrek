#include "decomposition.h"

Decomposition::Decomposition(SnpIndex *snpIndex, std::vector<Snp*> *snpList, Linkage *linkageMatrix, size_t thread, Region *regionInfo):m_snpIndex(snpIndex), m_snpList(snpList), m_linkage(linkageMatrix), m_thread(thread), m_regionInfo(regionInfo){}

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
    //DecompositionThread::checking = Eigen::MatrixXd::Zero(processSize, processSize); //DEBUG
    Eigen::VectorXd chiSq = Eigen::VectorXd::Zero(processSize);
    Eigen::VectorXd betaEstimate = Eigen::VectorXd::Zero(processSize);
	for(size_t i=0;i < processSize; ++i){
        if(snpLoc[i] < 0){
			throw "ERROR! Undefined behaviour. The location index is negative!";
        }
        else{
			betaEstimate(i) = (*m_snpList)[snpLoc[i]]->Getbeta();
			chiSq(i) = (*m_snpList)[snpLoc[i]]->GetsignedSqrtChiSq();
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
		//std::cerr << "Single thread" << std::endl;
		Eigen::MatrixXd variance = Eigen::MatrixXd::Zero(processSize, processSize);
		Eigen::MatrixXd additionVariance = Eigen::MatrixXd::Zero(processSize, processSize);
		Eigen::VectorXd result;
		result = m_linkage->solve(0, processSize, &betaEstimate, &chiSq, &variance, &additionVariance,Snp::GetmaxSampleSize());
		size_t copyStart = 0;
		size_t endOfProcess = blockSize /3*2;
		if(!chromosomeStart){
            copyStart = blockSize/3;
		}
		if(endOfProcess > processSize){ //If we encountered this, it means that it is the absolute end.
            endOfProcess = processSize;
		}
        for(size_t i=copyStart; i < endOfProcess; ++i){ /** I changed processSize to endOfProcess here */
			(*m_snpList)[snpLoc[i]]->Setheritability(result(i));
			(*m_snpList)[snpLoc[i]]->Setvariance(variance(i,i)); //The diagonal of the matrix contains the per snp variance
			(*m_snpList)[snpLoc[i]]->SetadditionVariance(additionVariance(i,i)); //The diagonal of the matrix contains the per snp variance

			for(size_t j = 0; j < processSize; ++j){ //For this snp, go through all the snp partners
				double covariance = variance(i,j);
				double additionCovariance = additionVariance(i,j);
				//DecompositionThread::checking(i,j) = variance(i,j);
				//DecompositionThread::checking(j,i) = variance(j,i);
				for(size_t regionIndex = 0; regionIndex < m_regionInfo->GetnumRegion(); ++regionIndex){
					if((*m_snpList)[snpLoc[i]]->GetFlag(regionIndex)&&
                        (*m_snpList)[snpLoc[j]]->GetFlag(regionIndex)){
                        m_regionInfo->Addvariance(covariance, regionIndex);
                        m_regionInfo->AddadditionVariance(additionCovariance, regionIndex);
					}
				}
			}
        }
        std::vector<double> regionBufferVariance(m_regionInfo->GetnumRegion(), 0.0);
        std::vector<double> regionBufferAdditionVariance(m_regionInfo->GetnumRegion(), 0.0);

        for(size_t i = endOfProcess; i < processSize; ++i){
            /** Doesn't matter with these three as if we are indeed not the last
             *  block of the chromosome, they will be re-wrote
             */
            (*m_snpList)[snpLoc[i]]->Setheritability(result(i));
			(*m_snpList)[snpLoc[i]]->Setvariance(variance(i,i));
			(*m_snpList)[snpLoc[i]]->SetadditionVariance(additionVariance(i,i));
            for(size_t j =  0; j < processSize; ++j){
                double covariance = variance(i,j);
                double additionCovariance = additionVariance(i,j);
				//DecompositionThread::checking(i,j) = variance(i,j);
				//DecompositionThread::checking(j,i) = variance(j,i);
                for(size_t regionIndex = 0; regionIndex < m_regionInfo->GetnumRegion(); ++regionIndex){
                    if((*m_snpList)[snpLoc[i]]->GetFlag(regionIndex)&&
                       (*m_snpList)[snpLoc[j]]->GetFlag(regionIndex)){
                        regionBufferVariance[regionIndex]+= covariance;
                        regionBufferAdditionVariance[regionIndex]+= additionCovariance;
                    }
                }
            }
        }
        for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->SetbufferVariance(regionBufferVariance[i],i);
            m_regionInfo->SetbufferAdditionVariance(regionBufferAdditionVariance[i], i);
        }

	}
	else if(stepSize > 0){
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
				if(i == startLoc.size()-1) garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate, &chiSq, m_linkage, &snpLoc, m_snpList, chromosomeStart, true, m_regionInfo));
				else garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate, &chiSq, m_linkage, &snpLoc, m_snpList, chromosomeStart, false, m_regionInfo));
				pthread_t *thread1 = new pthread_t();
				threadList.push_back(thread1);
				int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection.back());
				if(threadStatus != 0){
					throw "Failed to spawn thread with status: "+std::to_string(threadStatus);
				}
            }
			for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
			for(size_t i = 0; i < threadList.size(); ++i) delete threadList[i];
            for(size_t i = 0; i < garbageCollection.size(); ++i) delete garbageCollection[i];
			garbageCollection.clear();
            threadList.clear();
            //Now we do the extra job
			for(size_t i=m_thread; i < startLoc.size(); ++i){
				if(i == startLoc.size()-1) garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate, &chiSq, m_linkage, &snpLoc, m_snpList, chromosomeStart, true, m_regionInfo));
				else garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate, &chiSq, m_linkage, &snpLoc, m_snpList, chromosomeStart, false, m_regionInfo));
				pthread_t *thread1 = new pthread_t();
				threadList.push_back(thread1);
				int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection.back());
				if(threadStatus != 0){
					throw "Failed to spawn thread with status: "+std::to_string(threadStatus);
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
				if(i == startLoc.size()-1) garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate,&chiSq,m_linkage, &snpLoc, m_snpList, chromosomeStart, true, m_regionInfo));
				else garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate, &chiSq, m_linkage, &snpLoc, m_snpList, chromosomeStart, false, m_regionInfo));
				pthread_t *thread1 = new pthread_t();
				threadList.push_back(thread1);
				int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection.back());
				if(threadStatus != 0){
					throw "Failed to spawn thread with status: "+std::to_string(threadStatus);
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
        throw "Undefined behaviour! Step size should never be negative as block size is positive";
	}

	return completed;
}

