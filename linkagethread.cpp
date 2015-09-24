#include "linkagethread.h"
std::mutex LinkageThread::mtx;


LinkageThread::LinkageThread(bool correction, size_t vStart, size_t vEnd, size_t hEnd,  boost::ptr_deque<Genotype> *genotype, std::deque<size_t> *ldLoc,Eigen::MatrixXd *ldMatrix):m_correction(correction), m_vStart(vStart), m_vEnd(vEnd), m_hEnd(hEnd), m_genotype(genotype), m_ldLoc(ldLoc), m_ldMatrix(ldMatrix){}

LinkageThread::~LinkageThread(){}

void LinkageThread::ldConstruction(){
    //Now need to start the work
    for(size_t i = m_vStart; i <= m_vEnd; ++i){ //As vEnd is inclusive, it is not thread safe, need the mutex it if it is the vEnd level
        if(m_vEnd==i) mtx.lock();
        (*m_ldMatrix)(i,i) = 1.0;
        for(size_t j = i+1; j< m_hEnd; ++j){
            //Now perform the processing
            double rSquare=0.0;
            double r =0.0;
            //if(fabs((int)((*m_snpList)[(*m_snpLoc)[i]]->Getbp()-(*m_snpList)[(*m_snpLoc)[j]]->Getbp())) <= 2000000){
            (*m_genotype)[i].GetbothR(&((*m_genotype)[j]),m_correction, r, rSquare);
            // }
            (*m_ldMatrix)(i,j) = rSquare;
            (*m_ldMatrix)(j,i) = rSquare;

        }
        if(m_vEnd==i) mtx.unlock();
    }
}


void *LinkageThread::buildLd(void *in){
    ((LinkageThread *) in)->ldConstruction();
    return nullptr;
}



