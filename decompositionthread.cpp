#include "decompositionthread.h"

std::mutex DecompositionThread::decomposeMtx;

DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, bool lastOfBlock):m_start(start), m_length(length), m_betaEstimate(betaEstimate), m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList), m_chrStart(chrStart), m_lastOfBlock(lastOfBlock){}

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
	Eigen::VectorXd effective = Eigen::VectorXd::Constant(m_length, 1.0);
	Eigen::VectorXd result = m_linkage->solve(m_start, m_length, m_betaEstimate, &effective);
	size_t copyStart = m_length/3;
	size_t copyEnd = m_length/3;
	if(m_chrStart && m_start == 0){
		copyStart = 0;
        copyEnd += m_length/3;
	}
	if(m_lastOfBlock && m_start+m_length >= m_snpLoc->size()) copyEnd += m_length/3;

	DecompositionThread::decomposeMtx.lock();
	double effectiveNumber = 0.0;
	for(size_t i = copyStart; i < copyStart+copyEnd; ++i){
		(*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability(result(i));
		effectiveNumber+=effective(i);
	}
	m_linkage->Seteffective(effectiveNumber);
	DecompositionThread::decomposeMtx.unlock();
}
