#include "linkage.h"

Linkage::Linkage(size_t thread):m_thread(thread){}
Linkage::~Linkage(){}
std::mutex Linkage::mtx;


void Linkage::Initialize(boost::ptr_deque<Genotype> &genotype, const size_t &prevResiduals, const boost::ptr_vector<Interval> &blockSize){
	if(genotype.empty()){
        std::runtime_error("Cannot build LD without genotypes");
	}
	//Use previous informations
    if(prevResiduals == 0){
        m_linkage = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
        //m_linkageSqrt = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
    }
    else{
		Eigen::MatrixXd temp = m_linkage.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
        //Currently we only need the R2 matrix, so we will ignore the R matrix
		//temp = m_linkageSqrt.bottomRightCorner(prevResiduals,prevResiduals);
		//m_linkageSqrt= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		//m_linkageSqrt.topLeftCorner(prevResiduals, prevResiduals) = temp;

    }
}

void Linkage::buildLd(bool correction, size_t vStart, size_t vEnd, size_t hEnd, boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &ldLoc){
    for(size_t i = vStart; i <= vEnd; ++i){ //As vEnd is inclusive, it is not thread safe, need the mutex it if it is the vEnd level
        if(vEnd==i || vStart==i) mtx.lock();
        m_linkage(i,i) = 1.0;
        for(size_t j = i+1; j< hEnd; ++j){
            //Now perform the processing
            double rSquare=0.0;
            double r =0.0;
            //if(fabs((int)((*m_snpList)[(*m_snpLoc)[i]]->Getbp()-(*m_snpList)[(*m_snpLoc)[j]]->Getbp())) <= 2000000){
            genotype[i].GetbothR(&(genotype[j]),correction, r, rSquare);
            // }
            m_linkage(i,j) = rSquare;
            m_linkage(j,i) = rSquare;

        }
        if(vEnd==i || vStart==i) mtx.unlock();
    }
}

void Linkage::Construct(boost::ptr_deque<Genotype> &genotype, const size_t &genotypeIndex, const size_t& remainedLD, const boost::ptr_vector<Interval> &blockInfo, bool correction, std::deque<size_t> &ldLoc){
	if(genotype.empty())    throw std::runtime_error("There is no genotype to work on");
    size_t startRange =  genotypeIndex;
    size_t endRange=0;
    size_t i = genotypeIndex;
    size_t range = m_thread;
    std::string currentChr = blockInfo[genotypeIndex].getChr();
    if(remainedLD==0)range +=2;
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
    std::vector<std::thread> threadStore;
    //Launch a group of threads
    size_t workCount = 0;
    for(size_t i = startRange; i <= endRange; ++i){
        if(i+2<=endRange){
            //Long sausage
            threadStore.push_back(std::thread(&Linkage::buildLd, this, correction, blockInfo[i].getStart(), blockInfo[i].getEnd(),blockInfo[i+2].getEnd(), std::ref(genotype), std::ref(ldLoc)));
        }
        else if(i+1 <= endRange){
            //Short sausage
            threadStore.push_back(std::thread(&Linkage::buildLd, this, correction, blockInfo[i].getStart(), blockInfo[i].getEnd(),blockInfo[i+1].getEnd(), std::ref(genotype), std::ref(ldLoc)));
        }
        else{
            //triangle
            threadStore.push_back(std::thread(&Linkage::buildLd, this, correction, blockInfo[i].getStart(), blockInfo[i].getEnd(),blockInfo[i].getEnd(), std::ref(genotype), std::ref(ldLoc)));
        }
        workCount++;
        if(workCount >=m_thread){
            break;
        }
    }
    //Join the threads first
    for (size_t i = 0; i < threadStore.size(); ++i) {
        threadStore[i].join();
    }
    threadStore.clear();
    for(size_t i = startRange+workCount; i <=endRange; ++i){
        if(i+2<=endRange){
            //Long sausage
            threadStore.push_back(std::thread(&Linkage::buildLd, this, correction, blockInfo[i].getStart(), blockInfo[i].getEnd(),blockInfo[i+2].getEnd(), std::ref(genotype), std::ref(ldLoc)));
        }
        else if(i+1 <= endRange){
            //Short sausage
            threadStore.push_back(std::thread(&Linkage::buildLd, this, correction, blockInfo[i].getStart(), blockInfo[i].getEnd(),blockInfo[i+1].getEnd(), std::ref(genotype), std::ref(ldLoc)));
        }
        else{
            //triangle
            threadStore.push_back(std::thread(&Linkage::buildLd, this, correction, blockInfo[i].getStart(), blockInfo[i].getEnd(),blockInfo[i].getEnd(), std::ref(genotype), std::ref(ldLoc)));
        }
    }
    for (size_t i = 0; i < threadStore.size(); ++i) {
        threadStore[i].join();
    }
    threadStore.clear();
}

void Linkage::print(){
    std::cout << m_linkage << std::endl;
}
