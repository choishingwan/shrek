#include "decompositionthread.h"

std::mutex DecompositionThread::decomposeMtx;
Eigen::MatrixXd DecompositionThread::checking;


DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::VectorXd const * const betaEstimate, Eigen::VectorXd const * const sqrtChiSq, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, bool lastOfBlock, Region *regionInfo):m_start(start), m_length(length), m_betaEstimate(betaEstimate), m_sqrtChiSq(sqrtChiSq),m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList), m_chrStart(chrStart), m_lastOfBlock(lastOfBlock), m_regionInfo(regionInfo){}

DecompositionThread::~DecompositionThread()
{}


void *DecompositionThread::ThreadProcesser(void *in){
    ((DecompositionThread *) in)->solve();
	return nullptr;
}
void DecompositionThread::fullProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
//This is the case where we have the full block to decompose at once
    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    for(size_t i = 0; i < m_length; ++i){
        //Don't bother, just take everything
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        for(size_t j = 0; j < (size_t)(*variance).cols(); ++j){
            double covariance = (*variance)(i,j);
            Linkage::m_testing((*m_snpLoc)[m_start+i], j) =covariance;
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
	decomposeMtx.lock();
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
        for(size_t j = 0; j < (size_t)(*variance).cols() && j < m_length; ++j){
            double covariance = (*variance)(i,j);
            Linkage::m_testing((*m_snpLoc)[m_start+i], j) =covariance;
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Now the extra bit
    /*
    for(size_t i = m_length/3*2; i < m_length; ++i){
        for(size_t j = i-(m_length/3*2); j < m_length/3; ++j){
            double covariance=(*variance)(i,j);
            Linkage::m_testing((*m_snpLoc)[m_start+i], (*m_snpLoc)[m_start+j]) =covariance;
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    */
	decomposeMtx.lock();
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->Addvariance(regionVariance.at(i), i);
     	}
    decomposeMtx.unlock();
}

void DecompositionThread::normalProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
    //whether if it is the second last block

	std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    for(size_t i =m_length/3; i < m_length/3*2; ++i){

        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        for(size_t j = 0; j < (size_t)(*variance).cols(); ++j){
            double covariance = (*variance)(i,j);
            Linkage::m_testing((*m_snpLoc)[m_start+i], j) =covariance;
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Now the extra bit
    /*
    for(size_t i = m_length/3*2; i < m_length; ++i){
        for(size_t j = i-(m_length/3*2); j < m_length/3; ++j){
            double covariance=2*(*variance)(i,j);
            Linkage::m_testing((*m_snpLoc)[m_start+i], (*m_snpLoc)[m_start+j]) =covariance;
            Linkage::m_testing((*m_snpLoc)[m_start+j], (*m_snpLoc)[m_start+i]) =covariance;
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    */
	decomposeMtx.lock();
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->Addvariance(regionVariance.at(i), i);
     	}
	decomposeMtx.unlock();
}

void DecompositionThread::endBlockProcess(Eigen::MatrixXd const * const variance, Eigen::VectorXd const *const result){
    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    std::vector<double> regionBufferVariance(m_regionInfo->GetnumRegion(), 0.0);
    //The most complicated case (sort of);
    //The trick is, the variance matrix will always be with the correct dimension (else our methods has already failed)


    size_t actualProcessSize = (*variance).rows();
    for(size_t i =m_length/3; i < m_length/3*2; ++i){
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        for(size_t j = 0; j < (size_t)(*variance).cols(); ++j){
            double covariance = (*variance)(i,j);
            Linkage::m_testing((*m_snpLoc)[m_start+i], j) =covariance;
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    //Only add the variance of the top right missing part
    /*
    for(size_t i = m_length/3*2; i < actualProcessSize; ++i){
        for(size_t j = i-(m_length/3*2); j < m_length/3; ++j){
            double covariance = 2*(*variance)(i,j);
            Linkage::m_testing((*m_snpLoc)[m_start+i], (*m_snpLoc)[m_start+j]) =covariance;
            Linkage::m_testing((*m_snpLoc)[m_start+j], (*m_snpLoc)[m_start+i]) =covariance;
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
    }
    */
    //Finally, the bit where we want to put in the buffer
    for(size_t i = m_length/3*2; i < actualProcessSize; ++i){
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        for(size_t j = 0; j < (size_t)(*variance).cols(); ++j){
            double covariance = (*variance)(i,j);
            Linkage::m_testing((*m_snpLoc)[m_start+i], j) =covariance;
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionBufferVariance[regionIndex]+= covariance;
            }

        }
    }

	decomposeMtx.lock();
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->SetlastVariance(regionVariance.at(i),i);
            m_regionInfo->SetbufferVariance(regionBufferVariance.at(i),i);
     	}
    decomposeMtx.unlock();
}

void DecompositionThread::solve(){
    //Here is where I need to change stuff
    size_t betaLength = (*m_betaEstimate).rows();
    size_t processLength = m_length;
    if(m_lastOfBlock) processLength=betaLength-m_start;
    //Eigen::MatrixXd variance = Eigen::MatrixXd::Zero(processLength, processLength);
    Eigen::MatrixXd variance;
	Eigen::VectorXd result;
    result = m_linkage->solve(m_start, processLength, m_betaEstimate, m_sqrtChiSq, &variance, Snp::GetmaxSampleSize(),(*m_snpLoc)[m_start]);

    if(m_chrStart && m_lastOfBlock) fullProcess(&variance, &result);
    else if(m_chrStart&& m_start==0) chromosomeStartProcess(&variance, &result);
    else if(m_lastOfBlock) endBlockProcess(&variance, &result);
    else normalProcess(&variance, &result);

}
