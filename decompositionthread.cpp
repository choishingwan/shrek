#include "decompositionthread.h"

std::mutex DecompositionThread::decomposeMtx;
Eigen::MatrixXd DecompositionThread::checking;


DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::VectorXd const * const betaEstimate, Eigen::VectorXd const * const sqrtChiSq, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, Region *regionInfo):m_start(start), m_length(length), m_betaEstimate(betaEstimate), m_sqrtChiSq(sqrtChiSq),m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList), m_chrStart(chrStart), m_regionInfo(regionInfo){}

DecompositionThread::~DecompositionThread()
{}


void *DecompositionThread::ThreadProcesser(void *in){
    ((DecompositionThread *) in)->solve();
	return nullptr;
}
void DecompositionThread::fullProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
//This is the case where we have the full block to decompose at once
    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    size_t actualSize = (*variance).cols();
    for(size_t i = 0; i < actualSize; ++i){
        //Don't bother, just take everything
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        for(size_t j = 0; j < actualSize; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
	decomposeMtx.lock();
	//std::cerr << "Full process" << std::endl;
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->Addvariance(regionVariance.at(i), i);
     	}
    decomposeMtx.unlock();
}


void DecompositionThread::chromosomeStartProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
    //we need to include the front

	std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    for(size_t i =0; i < m_length/3*2; ++i){
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
    }
    //middle part
    for(size_t i = 0; i < m_length; ++i){
        for(size_t j = 0; j < m_length; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Top part
    for(size_t i = 0; i < m_length/3; ++i){
        for(size_t j = 0; j < m_length/3*2+i; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Bottom part
    for(size_t i = m_length/3*2; i < m_length; ++i){
        for(size_t j = i-m_length/3*2; j < m_length/3; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
	decomposeMtx.lock();
	//std::cerr << "Start of block" << std::endl;
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->Addvariance(regionVariance.at(i), i);
     	}
    decomposeMtx.unlock();
}

void DecompositionThread::normalProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
    //whether if it is the second last block

	std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
	//Middle part
    for(size_t i =m_length/3; i < m_length/3*2; ++i){
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        for(size_t j = 0; j < m_length; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Now the triangles
    for(size_t i = 0; i < m_length/3; ++i){
        for(size_t j = m_length/3*2; j < m_length/3*2+i+1; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    for(size_t i = m_length/3*2; i < m_length; ++i){
        for(size_t j = i-m_length/3*2; j < m_length/3; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
	decomposeMtx.lock();
	//std::cerr << "Normal parts" << std::endl;
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->Addvariance(regionVariance.at(i), i);
     	}
	decomposeMtx.unlock();
}

void DecompositionThread::endBlockProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);

    size_t actualProcessSize = (*variance).rows();
    for(size_t i =m_length/3; i < m_length/3*2; ++i){
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
    }
    //Middle part
    for(size_t i=m_length/3; i < m_length/3*2; ++i){
        for(size_t j = 0; j< actualProcessSize; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //bottom part
    for(size_t i = m_length/3*2; i < actualProcessSize; ++i){
        for(size_t j = i-m_length/3*2; j < actualProcessSize; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Top part
    for(size_t i = 0; i < m_length/3; ++i){
        for(size_t j = m_length/3*2; j < m_length/3*2+i+1; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }

	decomposeMtx.lock();
	//std::cerr << "End block"<< std::endl;
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->Addvariance(regionVariance.at(i),i);
     	}
    decomposeMtx.unlock();
}

void DecompositionThread::solve(){
    //Here is where I need to change stuff
    size_t betaLength = (*m_betaEstimate).rows();
    size_t processLength = m_length;
    if(betaLength-m_start-m_length < m_length/3) processLength = betaLength-m_start;
    Eigen::MatrixXd variance;
	Eigen::VectorXd result;
	decomposeMtx.lock();
    result = m_linkage->solve(m_start, processLength, m_betaEstimate, m_sqrtChiSq, &variance, Snp::GetmaxSampleSize(),(*m_snpLoc)[m_start]);
    decomposeMtx.unlock();
    size_t first = (*m_snpList)[(*m_snpLoc)[m_start]]->GetblockInfo();
    size_t last = (*m_snpList)[(*m_snpLoc)[m_start+processLength-1]]->GetblockInfo();
    if(first == 1 && last == 1) fullProcess(&variance, &result);
    else if(first==1 && last != 1) chromosomeStartProcess(&variance, &result);
    else if(last == 1 && first != 1) endBlockProcess(&variance, &result);
    else normalProcess(&variance, &result);
/*
    if(m_chrStart && m_lastOfBlock) fullProcess(&variance, &result);
    else if(m_chrStart&& m_start==0) chromosomeStartProcess(&variance, &result);
    else if(m_lastOfBlock) endBlockProcess(&variance, &result);
    else normalProcess(&variance, &result);
*/
}
