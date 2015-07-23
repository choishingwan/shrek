#include "decompositionthread.h"

std::mutex DecompositionThread::decomposeMtx;
Eigen::MatrixXd DecompositionThread::checking;


DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::VectorXd const * const betaEstimate, Eigen::VectorXd const * const sqrtChiSq, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, bool lastOfBlock, bool secondLastOfBlock,  Region *regionInfo):m_start(start), m_length(length), m_betaEstimate(betaEstimate), m_sqrtChiSq(sqrtChiSq),m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList), m_chrStart(chrStart), m_lastOfBlock(lastOfBlock), m_secondLastOfBlock(secondLastOfBlock), m_regionInfo(regionInfo){}

DecompositionThread::~DecompositionThread()
{}


void *DecompositionThread::ThreadProcesser(void *in){
    ((DecompositionThread *) in)->solve();
	return nullptr;
}

void DecompositionThread::chromosomeStartProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
    //we need to include the front

	decomposeMtx.lock();
    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    for(size_t i =0; i < m_length/3*2; ++i){
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        for(size_t j = 0; j < m_length; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Now the extra bit
    for(size_t i = m_length/3*2; i < m_length; ++i){
        for(size_t j = 0; j < m_length/3; ++j){
            double covariance=(*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
	//decomposeMtx.lock();
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->Addvariance(regionVariance.at(i), i);
     	}
    std::cerr << "Start of chromosome " << "\t" << (*m_snpLoc)[m_start] << "\t" << regionVariance.at(0)<< std::endl;
	decomposeMtx.unlock();
}

void DecompositionThread::normalProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
    //whether if it is the second last block

	decomposeMtx.lock();
    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    if((*m_snpLoc)[m_start]==596){
        std::cout << m_linkage->m_linkage.block(m_start, m_start, m_length, m_length) << std::endl;
    }
    for(size_t i =m_length/3; i < m_length/3*2; ++i){
        if((*m_snpLoc)[m_start]==596){
            std::cerr << (*m_snpList)[(*m_snpLoc)[m_start+i]]->GetrsId() << "\t" <<(*m_snpList)[(*m_snpLoc)[m_start+i]]->Getheritability() << "\t" << (*result)(i) << std::endl;

        }
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        for(size_t j = 0; j < m_length; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Now the extra bit
    for(size_t i = m_length/3*2; i < m_length; ++i){
        for(size_t j = 0; j < m_length/3; ++j){
            double covariance=2*(*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
	//decomposeMtx.lock();
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            if(m_secondLastOfBlock){
                m_regionInfo->AddbufferVariance(i,regionVariance.at(i));
            }
            else{
                m_regionInfo->Addvariance(regionVariance.at(i), i);
            }
     	}

	std::cerr << "Normal process " << (*m_snpLoc)[m_start] << "\t" << regionVariance.at(0) << "\t" << m_secondLastOfBlock << std::endl;
	decomposeMtx.unlock();
}

void DecompositionThread::endBlockProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    //The most complicated case (sort of);
    //The trick is, the variance matrix will always be with the correct dimension (else our methods has already failed)

	decomposeMtx.lock();

    size_t actualProcessSize = (*variance).rows();
    for(size_t i =m_length/3; i < actualProcessSize; ++i){
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        for(size_t j = 0; j < actualProcessSize; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Only add the variance of the top right missing part
    for(size_t i = 0; i < m_length/3; ++i){
        for(size_t j = m_length/3*2; j < actualProcessSize; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
	//decomposeMtx.lock();
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->AddbufferVariance(i,regionVariance.at(i));
     	}
    std::cerr << "Last of block" << "\t" << (*m_snpLoc)[m_start] << "\t" << actualProcessSize << "\t" << regionVariance.at(0)<< std::endl;
	decomposeMtx.unlock();
}

void DecompositionThread::solve(){
    //Here is where I need to change stuff
    size_t betaLength = (*m_betaEstimate).rows();
    size_t processLength = m_length;
    if(m_lastOfBlock) processLength=betaLength-m_start;
    Eigen::MatrixXd variance = Eigen::MatrixXd::Zero(processLength, processLength);
	Eigen::VectorXd result;
    result = m_linkage->solve(m_start, processLength, m_betaEstimate, m_sqrtChiSq, &variance, Snp::GetmaxSampleSize());;

    if(m_chrStart&& m_start==0) chromosomeStartProcess(&variance, &result);
    else if(m_lastOfBlock) endBlockProcess(&variance, &result);
    else normalProcess(&variance, &result);
}
