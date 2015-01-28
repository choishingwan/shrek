#include "linkage.h"

Linkage::Linkage(size_t thread):m_thread(thread)
{
	//ctor
	m_numItem= 0;
    m_effectiveNumber=0.0;
}

Linkage::~Linkage()
{
	//dtor
}


void Linkage::triangularThread( const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype){
    //Make the thread region
    std::vector<LinkageThread*> garbageCollection;
    for(size_t i = 0; i < m_thread; ++i){
        garbageCollection.push_back(new LinkageThread(correction, endBlock, &m_effectiveNumber, &m_linkage, &genotype));
    }
    for(size_t i = startBlock; i < endBlock; ++i){
        garbageCollection[i%m_thread]->Addstart(i); //So we distribute the items to the corresponding "thread"
    }
    std::vector<pthread_t*> threadList;
    for(size_t i = 0; i < m_thread; ++i){
        pthread_t *thread1 = new pthread_t();
        threadList.push_back(thread1);
        int threadStatus = pthread_create( thread1, NULL, &LinkageThread::triangularProcess, garbageCollection[i]);
        if(threadStatus != 0){
            std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
            exit(-1);
        }
    }
    for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter){
		pthread_join(*threadList[threadIter], NULL);
    }

    for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter){
        delete(threadList[threadIter]);
	}
	threadList.clear();
    for(size_t cleaning = 0; cleaning < garbageCollection.size(); ++cleaning){
        delete(garbageCollection[cleaning]);
    }
    garbageCollection.clear();
}


void Linkage::rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype){
    std::vector<LinkageThread*> garbageCollection;
    std::vector<pthread_t*> threadList;
    size_t snpStart = start-height;
    size_t snpEnd = start;
    if(m_thread > (snpEnd-snpStart)){
        //if there are more threads than items, each will do one horizontal row
        for(size_t i = snpStart; i < snpEnd; ++i){
            garbageCollection.push_back(new LinkageThread(correction, i, i+1, start, width, &m_effectiveNumber, &m_linkage, &genotype));
            pthread_t *thread1 = new pthread_t();
            threadList.push_back(thread1);
            int threadStatus = pthread_create( thread1, NULL, &LinkageThread::rectangularProcess, garbageCollection[i]);
            if(threadStatus != 0){
                std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
                exit(-1);
            }
        }
    }
    else{
        size_t step = (snpEnd -snpStart ) / m_thread;
        size_t remaining = (snpEnd -snpStart )%m_thread;
        size_t current = snpStart;
        for(size_t i = 0; i < m_thread; ++i){
            if(remaining > 0){
                garbageCollection.push_back(new LinkageThread(correction, current, current+step+1, start, width, &m_effectiveNumber, &m_linkage, &genotype));
                remaining--;
                current++;
            }
            else{
                garbageCollection.push_back(new LinkageThread(correction, current, current+step, start, width,&m_effectiveNumber, &m_linkage, &genotype));
            }
            pthread_t *thread1 = new pthread_t();
            threadList.push_back(thread1);
            int threadStatus = pthread_create( thread1, NULL, &LinkageThread::rectangularProcess, garbageCollection[i]);
            if(threadStatus != 0){
                std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
                exit(-1);
            }
            current += step;
        }
    }
    for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter){
		pthread_join(*threadList[threadIter], NULL);
    }

    for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter){
        delete(threadList[threadIter]);
	}
    for(size_t i = 0; i < garbageCollection.size(); ++i){
        delete garbageCollection[i];
    }
}


void Linkage::Construct(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize, bool correction){
	if(genotype.empty()){
        return;
	}
	if(blockSize == 0){
        //Doesn't have to build the LD matrix when the block size is 0
        return;
	}
	//Use previous informations
    if(prevResiduals == 0){
        m_linkage = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
    }
    else{
		Eigen::MatrixXd temp = m_linkage.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
    }
    //Construct the LD
	std::vector<int> startLoc;
	size_t currentBlockSize = blockSize;
	if(currentBlockSize > genotype.size()){
        currentBlockSize = genotype.size();
	}
    size_t stepSize = currentBlockSize/3;
    //If stepsize == 0 -> Block size < 3, just do the whole thing directly

    if(stepSize == 0){
        for(size_t i = 0; i < genotype.size(); ++i){
            m_linkage(i,i) = 1.0;
            m_effectiveNumber += 1.0;
            for(size_t j = i+1; j < genotype.size(); ++j){
                double rSquare = genotype[i]->Getr(genotype[j], correction);
                if(rSquare > 1.0){
                    std::cerr << "Warning: rSquare bigger than 1" << std::endl;
                    rSquare = 1.0;
                }
                if(rSquare < 0.0){
                    std::cerr << "Warning: rSquare smaller than 0" << std::endl;
                    rSquare = 0.0;
                }
                m_linkage(i,j) = rSquare;
                m_linkage(j,i) = rSquare;
                m_effectiveNumber += 2.0*rSquare;
            }
        }
        std::cerr << m_effectiveNumber << std::endl;
    }
    else{
        size_t counting= 0;
        for(size_t i = 0; i < genotype.size(); i+= stepSize){
            if(counting < 2){
                if(prevResiduals==0){
                    startLoc.push_back(i);
                }
            }
            else{
                startLoc.push_back(i);
            }
            counting++;
        }
        if(!startLoc.empty() && genotype.size() - startLoc.back() < stepSize) startLoc.pop_back(); //If the remaining is less than a stepSize from the end
        if((startLoc.empty() || startLoc.size()==0) && genotype.size()!= 0 ){
            //When this is the only block left
            triangularThread(stepSize*2, genotype.size(), correction, genotype);
            rectangularThread(0, genotype.size(),stepSize*2, correction, genotype);
        }
        else{
            for(size_t i = 0; i  < startLoc.size()-1; ++i){
                triangularThread(startLoc[i], startLoc[i+1], correction, genotype);
            }
            if(startLoc.size() != 0 && !startLoc.empty()){ //Double safe
                triangularThread(startLoc.back(), genotype.size(), correction, genotype);
            }
            for(size_t i = 0; i < startLoc.size()-1; ++i){
                if(prevResiduals==0 && i == 0){

                } //do nothing
                else if(prevResiduals==0 && i==1){
                    //The smaller chunk
                    rectangularThread(startLoc[i], startLoc[i+1], stepSize, correction, genotype);
                }
                else{
                    rectangularThread(startLoc[i], startLoc[i+1], stepSize*2, correction, genotype);
                }
            }
            if(startLoc.size() >= 1){
                rectangularThread(startLoc.back(),genotype.size() , stepSize*2, correction, genotype);
            }
        }
        startLoc.clear();
    }

}
