#include "decompositionthread.h"

std::mutex DecompositionThread::decomposeMtx;

DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::VectorXd const * const betaEstimate, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, bool lastOfBlock):m_start(start), m_length(length), m_betaEstimate(betaEstimate),m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList), m_chrStart(chrStart), m_lastOfBlock(lastOfBlock){}

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
	Eigen::VectorXd result = m_linkage->solveChi(m_start, processLength, m_betaEstimate, variance, Snp::m_maxSampleSize);

	size_t copyStart = m_length/3;
	size_t copyEnd = m_length/3;
	if(m_chrStart && m_start == 0){
		copyStart = 0;
        copyEnd += m_length/3;
	}
	if(m_lastOfBlock) copyEnd = processLength-copyStart;
    for(size_t i = copyStart; i < copyStart+copyEnd; ++i){
		(*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability(result(i));
		(*m_snpList)[(*m_snpLoc)[m_start+i]]->Setvariance(variance(i,i)); //The diagonal of the matrix contain the per snp variance
        for(size_t j = 0; j < processLength; ++j){ //For this snp, go through all the snp partners
			(*m_snpList)[(*m_snpLoc)[m_start+j]]->Setheritability(result(i));
        }


	}
}
