#include "decompositionthread.h"

std::mutex DecompositionThread::decomposeMtx;

DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::VectorXd const * const betaEstimate, Eigen::VectorXd const * const signValue, Eigen::VectorXd const * const chiSq, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, bool lastOfBlock):m_start(start), m_length(length), m_betaEstimate(betaEstimate), m_signValue(signValue),m_chiSq(chiSq),m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList), m_chrStart(chrStart), m_lastOfBlock(lastOfBlock){}

DecompositionThread::~DecompositionThread()
{}


void *DecompositionThread::ThreadProcesser(void *in){
    ((DecompositionThread *) in)->solve();
	return nullptr;
}

void DecompositionThread::solve(){
    size_t betaLength = (*m_betaEstimate).rows();
    size_t processLength = m_length;
    if(m_lastOfBlock) processLength=betaLength-m_start;
    Eigen::MatrixXd variance = Eigen::MatrixXd::Zero(processLength, processLength);
    Eigen::MatrixXd additionVariance = Eigen::MatrixXd::Zero(processLength, processLength);
	Eigen::VectorXd result = m_linkage->solveChi(m_start, processLength, m_betaEstimate, m_signValue, m_chiSq, &variance, &additionVariance, Snp::m_maxSampleSize);

	size_t copyStart = m_length/3;
	size_t copyEnd = m_length/3;
	if(m_chrStart && m_start == 0){
		copyStart = 0;
        copyEnd += m_length/3;
	}
	if(m_lastOfBlock) copyEnd = processLength-copyStart;
    std::vector<double> regionVariance(Region::regionVariance.size(), 0.0);
    std::vector<double> regionAdditionVariance(Region::regionAdditionVariance.size(), 0.0);
    for(size_t i = copyStart; i < copyStart+copyEnd; ++i){
		(*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability(result(i));
		(*m_snpList)[(*m_snpLoc)[m_start+i]]->Setvariance(variance(i,i)); //The diagonal of the matrix contain the per snp variance
		(*m_snpList)[(*m_snpLoc)[m_start+i]]->SetadditionVariance(additionVariance(i,i)); //The diagonal of the matrix contain the per snp variance

        for(size_t j = 0; j < processLength; ++j){ //For this snp, go through all the snp partners
			double covariance = variance(i,j);
			double additionCovariance = additionVariance(i,j);
			for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                if((*m_snpList)[(*m_snpLoc)[m_start+i]]->GetFlag(regionIndex)&&(*m_snpList)[(*m_snpLoc)[m_start+j]]->GetFlag(regionIndex)){
					regionVariance[regionIndex]+= covariance;
					regionAdditionVariance[regionIndex]+= additionCovariance;
                }
			}
        }
	}
	decomposeMtx.lock();
		for(size_t i = 0; i < Region::regionVariance.size(); ++i){
			Region::regionVariance[i]+=regionVariance[i];
			Region::regionAdditionVariance[i]+=regionAdditionVariance[i];
		}
	decomposeMtx.unlock();

}
