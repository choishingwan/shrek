#include "decomposition.h"

//require major revision. There are BUGS here. So we need to make it robust to changes.

Decomposition::Decomposition( std::vector<Snp*> *snpList, Linkage *linkageMatrix, size_t thread, Region *regionInfo):m_snpList(snpList), m_linkage(linkageMatrix), m_thread(thread), m_regionInfo(regionInfo){}

Decomposition::~Decomposition()
{
	//dtor
}

ProcessCode Decomposition::Decompose(const size_t &blockSize, std::deque<size_t> &snpLoc, std::deque<Genotype*> &genotype, bool chromosomeStart, bool chromosomeEnd){
    //First, build the vector for beta
    if(genotype.size() == 0){
        //new algorithm should not have any situation where genotype size is 0
        throw "No genotype to work on";
    }
    //processSize by default should be the number of Snps
    size_t processSize = snpLoc.size();
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
        //Just in case if the blockSize is much bigger than the number of Snps
        currentBlockSize = processSize;
	}
    std::vector<size_t> startLoc;
	size_t stepSize = currentBlockSize/3;
    //Now we need to obtain the actual work coordinates
    //For the last two block, their variance should always put into buffer
    //Other than that, there is no changes.
    //Also, ignore the single thread situation to reduce the size of the code



    if(stepSize == 0){
        //only possible if the blockSize is less than 3 e.g. 2 Snps or 1 Snps
        //The only time when this is possible is when there is only 2 or 1 Snp(s)
        //in the input file.
        //Will still treat it specially.
        if(processSize > 2){
            std::cerr << "Something is wrong. I don't expect step size = 0 when there are more than 2 snp input" << std::endl;
            throw "Unexpected error";
        }
        //Eigen::MatrixXd variance = Eigen::MatrixXd::Zero(processSize, processSize);
        Eigen::MatrixXd variance;
		Eigen::VectorXd result;
		result = m_linkage->solve(0, processSize, &betaEstimate, &chiSq, &variance, Snp::GetmaxSampleSize(),snpLoc[0]);
        for(size_t i = 0; i < processSize; ++i){
            (*m_snpList).at(snpLoc.at(i))->Setheritability(result(i));
            (*m_snpList).at(snpLoc.at(i))->Setvariance(variance(i,i));
            for(size_t j = 0; j < processSize;++j){
                double covariance = variance(i,j);
                for(size_t regionIndex = 0; regionIndex < m_regionInfo->GetnumRegion(); ++regionIndex){
					if((*m_snpList)[snpLoc[i]]->GetFlag(regionIndex)&&
                        (*m_snpList)[snpLoc[j]]->GetFlag(regionIndex)){
                        m_regionInfo->Addvariance(covariance, regionIndex);
					}
				}
            }
        }

    }
    else if(stepSize> 0){
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
        for(size_t i = 0; i < startLoc.size(); ++i){
            garbageCollection.push_back(new DecompositionThread(startLoc[i], currentBlockSize, &betaEstimate,&chiSq,m_linkage, &snpLoc, m_snpList, chromosomeStart, m_regionInfo));
        }
        if(m_thread >= garbageCollection.size()){
            //More thread than work
            for(size_t i = 0; i < garbageCollection.size(); ++i){
                pthread_t *thread1 = new pthread_t();
				threadList.push_back(thread1);
				int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection[i]);
				if(threadStatus != 0){
					throw "Failed to spawn thread with status: "+std::to_string(threadStatus);
				}
            }
			for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
			for(size_t i = 0; i < threadList.size(); ++i) delete threadList[i];
            for(size_t i = 0; i < garbageCollection.size(); ++i) delete garbageCollection[i];
			garbageCollection.clear();
            threadList.clear();
        }
        else{
            //More work than thread
            for(size_t i = 0; i < m_thread; ++i){
                pthread_t *thread1 = new pthread_t();
				threadList.push_back(thread1);
				int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection[i]);
				if(threadStatus != 0){
					throw "Failed to spawn thread with status: "+std::to_string(threadStatus);
				}
            }
			for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
			for(size_t i = 0; i < threadList.size(); ++i) delete threadList[i];
			threadList.clear();
            for(size_t i= m_thread; i < garbageCollection.size(); ++i){
                pthread_t *thread1 = new pthread_t();
				threadList.push_back(thread1);
				int threadStatus = pthread_create( thread1, NULL, &DecompositionThread::ThreadProcesser, garbageCollection[i]);
				if(threadStatus != 0){
					throw "Failed to spawn thread with status: "+std::to_string(threadStatus);
				}
            }

			for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
			for(size_t i = 0; i < threadList.size(); ++i) delete threadList[i];
            for(size_t i = 0; i < garbageCollection.size(); ++i) delete garbageCollection[i];
			garbageCollection.clear();
            threadList.clear();

        }

    }
    else{
        throw "Undefined behaviour! Step size should never be negative as block size is positive";
	}

    return completed;
}

