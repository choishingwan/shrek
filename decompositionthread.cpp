#include "decompositionthread.h"

std::mutex DecompositionThread::decomposeMtx;
Eigen::MatrixXd DecompositionThread::checking;


DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::VectorXd const * const betaEstimate, Eigen::VectorXd const * const sqrtChiSq, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, Region *regionInfo):m_start(start), m_length(length), m_betaEstimate(betaEstimate), m_sqrtChiSq(sqrtChiSq),m_sampleMatrix(nullptr),m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList), m_chrStart(chrStart), m_regionInfo(regionInfo){}

DecompositionThread::~DecompositionThread()
{}

void *DecompositionThread::ThreadProcesser(void *in){
    ((DecompositionThread *) in)->solve();
	return nullptr;
}


void *DecompositionThread::SampleProcesser(void *in){
    ((DecompositionThread *) in)->sampleSolve();
    return nullptr;
}



void DecompositionThread::fullProcess(Eigen::VectorXd const * const perSnpEffect, Eigen::VectorXd const *const result, Eigen::VectorXd const *const effectiveReturnResult){
//This is the case where we have the full block to decompose at once
    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    size_t actualSize = (*result).rows();
    for(size_t i = 0; i < actualSize; ++i){
        //Don't bother, just take everything
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->SeteffectiveNumber((*effectiveReturnResult)(i));
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->SetsnpLDSC((*perSnpEffect)(i));
        /*
        for(size_t j = 0; j < actualSize; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
        */
    }
    //std::cerr << "Full process" << std::endl;
    /*
	decomposeMtx.lock();
	//std::cerr << "Full process" << std::endl;
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->Addvariance(regionVariance.at(i), i);
     	}
    decomposeMtx.unlock();
    */
}


void DecompositionThread::chromosomeStartProcess(Eigen::VectorXd const * const perSnpEffect, Eigen::VectorXd const *const result, Eigen::VectorXd const *const effectiveReturnResult){
    //we need to include the front

	std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    for(size_t i =0; i < m_length/3*2; ++i){
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->SeteffectiveNumber((*effectiveReturnResult)(i));
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->SetsnpLDSC((*perSnpEffect)(i));
    }
    //std::cerr << "Start of block" << std::endl;
    /*
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
    */
}

void DecompositionThread::normalProcess(Eigen::VectorXd const * const perSnpEffect, Eigen::VectorXd const *const result, Eigen::VectorXd const *const effectiveReturnResult){
    //whether if it is the second last block

	std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
	//Middle part
    for(size_t i =m_length/3; i < m_length/3*2; ++i){
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->SeteffectiveNumber((*effectiveReturnResult)(i));
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->SetsnpLDSC((*perSnpEffect)(i));

        /*
        for(size_t j = 0; j < m_length; ++j){
            double covariance = (*variance)(i,j);
            for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                regionVariance[regionIndex]+= covariance;
            }
        }
        */
    }
    //std::cerr << "Normal parts" << std::endl;
    /*
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
	*/
}

void DecompositionThread::endBlockProcess(Eigen::VectorXd const * const perSnpEffect, Eigen::VectorXd const *const result, Eigen::VectorXd const *const effectiveReturnResult){
    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);

    size_t actualProcessSize = (*result).rows();
    for(size_t i =m_length/3; i < actualProcessSize; ++i){
        //std::cerr << "Set: " << (*m_snpLoc)[m_start+i] << "\t" << (*result)(i) << std::endl;
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability((*result)(i));
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->SeteffectiveNumber((*effectiveReturnResult)(i));
        (*m_snpList)[(*m_snpLoc)[m_start+i]]->SetsnpLDSC((*perSnpEffect)(i));
    }
    //std::cerr << "End block"<< std::endl;
    /*
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
    */
}

void DecompositionThread::solve(){
    //Here is where I need to change stuff
    size_t betaLength = (*m_betaEstimate).rows();
    size_t processLength = m_length;
    if(betaLength-m_start-m_length < m_length/3) processLength = betaLength-m_start;
    Eigen::VectorXd perSnpEffect;
	Eigen::VectorXd result;
	Eigen::VectorXd  effectiveReturnResult;

    result = m_linkage->solve(m_start, processLength, m_betaEstimate, m_sqrtChiSq, &perSnpEffect, &effectiveReturnResult,
 Snp::GetmaxSampleSize(),(*m_snpLoc)[m_start]);

    size_t first = (*m_snpList)[(*m_snpLoc)[m_start]]->GetblockInfo();
    size_t last = (*m_snpList)[(*m_snpLoc)[m_start+processLength-1]]->GetblockInfo();
    if(first == 1 && last == 1) fullProcess(&perSnpEffect, &result,&effectiveReturnResult);
    else if(first==1 && last != 1) chromosomeStartProcess(&perSnpEffect, &result,&effectiveReturnResult);
    else if(last == 1 && first != 1) endBlockProcess(&perSnpEffect, &result,&effectiveReturnResult);
    else normalProcess(&perSnpEffect, &result,&effectiveReturnResult);
/*
    if(m_chrStart && m_lastOfBlock) fullProcess(&variance, &result);
    else if(m_chrStart&& m_start==0) chromosomeStartProcess(&variance, &result);
    else if(m_lastOfBlock) endBlockProcess(&variance, &result);
    else normalProcess(&variance, &result);
*/
}

/**
 *  This sections are for risk prediction
 */


DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::MatrixXd const * const sampleMatrix,Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, std::vector<double> *samplePheno, bool chrStart):m_start(start), m_length(length), m_betaEstimate(nullptr), m_sqrtChiSq(nullptr), m_sampleMatrix(sampleMatrix),m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList),m_samplePheno(samplePheno), m_chrStart(chrStart){}

void DecompositionThread::sampleSolve(){
    size_t betaLength = (*m_sampleMatrix).rows();
    size_t processLength = m_length;
    if(betaLength-m_start-m_length < m_length/3) processLength = betaLength-m_start;
	Eigen::MatrixXd result;
    result = m_linkage->solve(m_start, processLength, m_sampleMatrix,Snp::GetmaxSampleSize(),(*m_snpLoc)[m_start]);
    size_t first = (*m_snpList)[(*m_snpLoc)[m_start]]->GetblockInfo();
    size_t last = (*m_snpList)[(*m_snpLoc)[m_start+processLength-1]]->GetblockInfo();
    if(first == 1 && last == 1) fullProcess(&result,m_samplePheno);
    else if(first==1 && last != 1) chromosomeStartProcess(&result,m_samplePheno);
    else if(last == 1 && first != 1) endBlockProcess(&result,m_samplePheno);
    else normalProcess(&result,m_samplePheno);

}

void DecompositionThread::chromosomeStartProcess(Eigen::MatrixXd const *const result, std::vector<double> *m_samplePheno){
    decomposeMtx.lock();
    for(size_t i = 0; i < (*m_samplePheno).size(); ++i){
        (*m_samplePheno)[i] += (*result).col(i).segment(0, m_length/3*2).sum();
    }
    decomposeMtx.unlock();
}

void DecompositionThread::endBlockProcess(Eigen::MatrixXd const *const result, std::vector<double> *m_samplePheno){
    size_t actualProcessSize = (*result).rows();
    decomposeMtx.lock();
    for(size_t i = 0; i < (*m_samplePheno).size(); ++i){
        (*m_samplePheno)[i] += (*result).col(i).segment(m_length/3,actualProcessSize-m_length/3).sum();
    }
    decomposeMtx.unlock();
}

void DecompositionThread::normalProcess(Eigen::MatrixXd const *const result, std::vector<double> *m_samplePheno){
    decomposeMtx.lock();
    for(size_t i = 0; i < (*m_samplePheno).size(); ++i){
        (*m_samplePheno)[i] += (*result).col(i).segment(m_length/3,m_length/3).sum();
    }
    decomposeMtx.unlock();
}

void DecompositionThread::fullProcess(Eigen::MatrixXd const *const result, std::vector<double> *m_samplePheno){
    decomposeMtx.lock();
    for(size_t i = 0; i < (*m_samplePheno).size(); ++i){
        //Don't bother, just take everything
        (*m_samplePheno)[i] += (*result).col(i).sum();
    }
    decomposeMtx.unlock();
}

