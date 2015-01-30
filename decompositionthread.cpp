#include "decompositionthread.h"

std::mutex DecompositionThread::decomposeMtx;

DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, bool chrEnd):m_start(start), m_length(length), m_betaEstimate(betaEstimate), m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList), m_chrStart(chrStart), m_chrEnd(chrEnd){}

DecompositionThread::~DecompositionThread()
{
	//dtor
}


void *DecompositionThread::ThreadProcesser(void *in){
    struct DecompositionThread *input = (DecompositionThread *) in;
    input->solve();
	return nullptr;
}

void DecompositionThread::solve(){
	DecompositionThread::decomposeMtx.lock();

    Eigen::VectorXd result = m_linkage->solve(m_start, m_length, m_betaEstimate);
	size_t copyStart = m_length/3;
	size_t copyEnd = m_length/3;
	if(m_chrStart && m_start == 0) copyStart = 0;
	if(m_chrEnd && m_start+m_length >= m_snpLoc->size()) copyEnd = m_length;


	for(size_t i = copyStart; i < copyEnd; ++i){
		(*m_snpList)[(*m_snpLoc)[i]]->Setheritability(result(i));
	}
	DecompositionThread::decomposeMtx.unlock();
}
