#include "linkagethread.h"

bool LinkageThread::Getcorrection() const { return m_correction; }
size_t LinkageThread::GetsnpStart() const { return m_snpStart; };
size_t LinkageThread::GetsnpEnd() const { return m_snpEnd; };
size_t LinkageThread::GetboundStart() const { return m_boundStart; }
size_t LinkageThread::GetboundEnd() const { return m_boundEnd; }
size_t LinkageThread::GetstartLoc(size_t i) const{ return m_startLoc[i]; }
size_t LinkageThread::GetsizeOfStart() const { return m_startLoc.size(); }
Eigen::MatrixXd *LinkageThread::Getld() { return m_ldMatrix; }
std::deque<Genotype* > *LinkageThread::Getgenotype(){ return m_genotype; }

void LinkageThread::Addstart(size_t i){
	m_startLoc.push_back(i);
}

LinkageThread::LinkageThread(bool correction, const size_t blockEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype):m_correction(correction), m_boundEnd(blockEnd), m_ldMatrix(ldMatrix), m_genotype(genotype){}

LinkageThread::LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype):m_correction(correction), m_snpStart(snpStart), m_snpEnd(snpEnd), m_boundStart(boundStart), m_boundEnd(boundEnd), m_ldMatrix(ldMatrix), m_genotype(genotype){}

LinkageThread::~LinkageThread(){}

void LinkageThread::triangularProcess(){
    size_t sizeOfItem = m_startLoc.size();
    for(size_t i = 0; i < sizeOfItem; ++i){
        size_t j = m_startLoc[i];
        (*m_ldMatrix)(j,j) = 1.0;
        for(size_t k = m_boundEnd-1; k > j; --k){
            if((*m_ldMatrix)(j,k) == 0.0){
                double rSquare=(*m_genotype)[j]->Getr( (*m_genotype)[k], m_correction);
                (*m_ldMatrix)(j, k) = rSquare;
                (*m_ldMatrix)(k,j) = rSquare;
            }
            else break;
        }
    }

}

void LinkageThread::rectangularProcess(){
    for(size_t i = m_snpStart; i < m_snpEnd; ++i){
        for(size_t j = m_boundEnd-1; j >= m_boundStart; ++j){
            if(j == i){
                (*m_ldMatrix)(i,i) = 1.0; //Let's just assume that it is duplicated
            }
            else if((*m_ldMatrix)(i,j) == 0.0){
                double rSquare = (*m_genotype)[i]->Getr((*m_genotype)[j],m_correction);
				(*m_ldMatrix)(i,j) = rSquare;
				(*m_ldMatrix)(j,i) = rSquare;
            }
            else break;
        }
    }


}


void *LinkageThread::triangularProcess(void *in){
    struct LinkageThread *input = (LinkageThread *) in;
    input->triangularProcess();
    return nullptr;
}

void *LinkageThread::rectangularProcess(void *in){
    struct LinkageThread *input = (LinkageThread *) in;
    input->rectangularProcess();
    return nullptr;
}


