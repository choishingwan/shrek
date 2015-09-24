#include "linkage.h"

Linkage::Linkage(size_t thread):m_thread(thread){}


void Linkage::Initialize(boost::ptr_deque<Genotype> &genotype, const size_t &prevResiduals, const boost::ptr_vector<Interval> &blockSize){
	if(genotype.empty()){
        std::runtime_error("Cannot build LD without genotypes")
	}
	//Use previous informations
    if(prevResiduals == 0){
        m_linkage = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
        m_linkageSqrt = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
    }
    else{
        blockSize
		Eigen::MatrixXd temp = m_linkage.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
        //Currently we only need the R2 matrix, so we will ignore the R matrix
		//temp = m_linkageSqrt.bottomRightCorner(prevResiduals,prevResiduals);
		//m_linkageSqrt= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		//m_linkageSqrt.topLeftCorner(prevResiduals, prevResiduals) = temp;

    }
}

void Linkage::buildLd(boost::ptr_deque<Genotype> &genotype,  size_t vStart, size_t hEnd, size_t vEnd, std::deque<size_t> &ldLoc){

}

void Linkage::Construct(boost::ptr_deque<Genotype> &genotype, const size_t &genotypeIndex, const size_t& remainedLD, const boost::ptr_vector<Interval> &blockSize, bool correction, std::deque<size_t> &ldLoc){
	if(genotype.empty())    throw std::runtime_error("There is no genotype to work on");
    size_t startRange =  genotypeIndex;
    size_t endRange=0;
    size_t i = genotypeIndex;
    size_t range = m_thread;
    if(remainedLd==0)range +=2;
    for(; i< genotypeIndex+range && i < blockInfo.size(); ++i){
            if(blockInfo[i].getChr().compare(currentChr)!=0){
            i= i-1; //This mean we working on the last block of this chromosome
            break;
        }
    }
    endRange = i;
    // So now the startRange and endRange will contain the index of the intervals to include in the LD construction
    // Each thread will process one sausage       \------------|
    //                                             \-----------|

    //Change this into threading
    std::vector<thread> threads;
    //Launch a group of threads
    size_t workCount = 0;
    for(size_t i = startRange; i <= endRange; ++i){
        if(i+2<=endRange){
            //Long sausage
            std::thread(buildLd, correction, blockSize[i].getStart(), blockSize[i].getEnd(),blockSize[i+2].getEnd(), genotype, ldLoc)
        }
        else if(i+1 <= endRange){
            //Short sausage
            std::thread(buildLd, correction, blockSize[i].getStart(), blockSize[i].getEnd(),blockSize[i+1].getEnd(), genotype, ldLoc)
        }
        else{
            //triangle
            std::thread(buildLd, correction, blockSize[i].getStart(), blockSize[i].getEnd(),blockSize[i].getEnd(), genotype, ldLoc)
        }
        workCount++;
        if(workCount >=m_thread){
            break;
        }
    }
    //Join the threads first
    for (size_t i = 0; i < threads.size(); ++i) {
        threads[i].join();
    }
    threads.clear();
    for(size_t i = startRange+workCount; i <=endRange; ++i){
        if(i+2<=endRange){
            //Long sausage
            std::thread(buildLd, correction, blockSize[i].getStart(), blockSize[i].getEnd(),blockSize[i+2].getEnd(), genotype, ldLoc)
        }
        else if(i+1 <= endRange){
            //Short sausage
            std::thread(buildLd, correction, blockSize[i].getStart(), blockSize[i].getEnd(),blockSize[i+1].getEnd(), genotype, ldLoc)
        }
        else{
            //triangle
            std::thread(buildLd, correction, blockSize[i].getStart(), blockSize[i].getEnd(),blockSize[i].getEnd(), genotype, ldLoc)
        }
    }
    for (size_t i = 0; i < threads.size(); ++i) {
        threads[i].join();
    }
    threads.clear();
    boost::ptr_vector<LinkageThread> garbageCollection;
    boost::ptr_vector<pthread_t> threadList;
    for(size_t i = startRange; i <= endRange; ++i){
        if(i+2<=endRange){
            //Long sausage
            garbageCollection.push_back(new LinkageThread(correction, blockSize[i].getStart(), blockSize[i].getEnd(),blockSize[i+2].getEnd(), &genotype, &ldLoc, &m_linkage));
        }
        else if(i+1 <= endRange){
            //Short sausage
            garbageCollection.push_back(new LinkageThread(correction, blockSize[i].getStart(), blockSize[i].getEnd(),blockSize[i+1].getEnd(), &genotype, &ldLoc, &m_linkage));
        }
        else{
            //triangle
            garbageCollection.push_back(new LinkageThread(correction, blockSize[i].getStart(), blockSize[i].getEnd(),blockSize[i].getEnd(), &genotype, &ldLoc, &m_linkage));
        }
    }
    //Now send jobs to threads
    if(m_thread >= garbageCollection.size()){
        //More thread than work
        for(size_t i = 0; i < garbageCollection.size(); ++i){
            pthread_t *thread1 = new pthread_t();
            threadList.push_back(thread1);
            int threadStatus = pthread_create( thread1, NULL, &LinkageThread::buildLd, garbageCollection[i]);
			if(threadStatus != 0){
				throw std::runtime_error("Failed to spawn thread with status: "+std::to_string(threadStatus));
            }
        }
        for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
        garbageCollection.clear();
        threadList.clear();
    }
    else{
        //More work than thread
        if(garbageCollection.size() > 2*m_thread){
            throw std::logic_error("Sorry, logic error again, somehow I will need more threads than I should");
        }
        for(size_t i = 0; i < m_thread; ++i){
            pthread_t *thread1 = new pthread_t();
            threadList.push_back(thread1);
            int threadStatus = pthread_create( thread1, NULL, &LinkageThread::buildLd, garbageCollection[i]);
            if(threadStatus != 0){
                throw std::runtime_error("Failed to spawn thread with status: "+std::to_string(threadStatus));
            }
        }
        for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
        threadList.clear();
        for(size_t i= m_thread; i < garbageCollection.size(); ++i){
            pthread_t *thread1 = new pthread_t();
            threadList.push_back(thread1);
            int threadStatus = pthread_create( thread1, NULL, &LinkageThread::buildLd, garbageCollection[i]);
            if(threadStatus != 0){
                throw std::runtime_error("Failed to spawn thread with status: "+std::to_string(threadStatus));
            }
        }

        for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter) pthread_join(*threadList[threadIter], NULL);
        garbageCollection.clear();
        threadList.clear();
    }
}
