#include "linkagethread.h"
std::mutex LinkageThread::mtx;

void LinkageThread::Addstart(size_t i){ m_startLoc.push_back(i); }

LinkageThread::LinkageThread(bool correction, const size_t blockEnd, Eigen::MatrixXd *ldMatrix, Eigen::MatrixXd *ldMatrixSqrt, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList):m_correction(correction), m_boundEnd(blockEnd), m_ldMatrix(ldMatrix), m_ldMatrixSqrt(ldMatrixSqrt),m_genotype(genotype), m_snpLoc(snpLoc), m_perfectLd(perfectLd), m_snpList(snpList){}

LinkageThread::LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, Eigen::MatrixXd *ldMatrix, Eigen::MatrixXd *ldMatrixSqrt, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList):m_correction(correction), m_snpStart(snpStart), m_snpEnd(snpEnd), m_boundStart(boundStart), m_boundEnd(boundEnd), m_ldMatrix(ldMatrix), m_ldMatrixSqrt(ldMatrixSqrt), m_genotype(genotype), m_snpLoc(snpLoc), m_perfectLd(perfectLd), m_snpList(snpList){}

LinkageThread::~LinkageThread(){}

void LinkageThread::triangularProcess(){
    /** m_boundEnd is the end, not inclusive */
    size_t numOfSnp = m_startLoc.size(); //Number of snp to work on
    std::vector<size_t> perfectLd;
    for(size_t i = 0; i < numOfSnp; ++i){
        size_t j = m_startLoc[i];
        (*m_ldMatrix)(j,j) = 1.0;
        (*m_ldMatrixSqrt)(j,j) = 1.0;

        for(size_t k = m_boundEnd-1; k > j; --k){
            if((*m_ldMatrix)(j,k) == 0.0){ //If 0.0, then it is new and has never been calculated before
                double rSquare= 0.0;
                double r = 0.0;
                (*m_genotype)[j]->GetbothR( (*m_genotype)[k], m_correction, r, rSquare);
                (*m_ldMatrix)(j, k) = rSquare;
                (*m_ldMatrix)(k,j) = rSquare;
                (*m_ldMatrixSqrt)(j, k) = r;
                (*m_ldMatrixSqrt)(k,j) = r;
            }
            else break;
        }
    }
    LinkageThread::mtx.lock();
    (*m_perfectLd).insert((*m_perfectLd).end(), perfectLd.begin(), perfectLd.end());
    LinkageThread::mtx.unlock();
}

void LinkageThread::rectangularProcess(){

    std::vector<size_t> perfectLd;
    for(size_t i = m_snpStart; i < m_snpEnd; ++i){
        for(size_t j = m_boundEnd-1; j >= m_boundStart; --j){
            if(j == i){
                (*m_ldMatrix)(i,i) = 1.0; //Let's just assume that it is duplicated
                (*m_ldMatrixSqrt)(i,i) = 1.0;
            }
            else if((*m_ldMatrix)(i,j) == 0.0){
                double rSquare=0.0;
                double r =0.0;
                (*m_genotype)[i]->GetbothR((*m_genotype)[j],m_correction, r, rSquare);
				(*m_ldMatrix)(i,j) = rSquare;
				(*m_ldMatrix)(j,i) = rSquare;
				(*m_ldMatrixSqrt)(i,j) = r;
				(*m_ldMatrixSqrt)(j,i) = r;

            }
            else break;
        }
    }
    LinkageThread::mtx.lock();
    (*m_perfectLd).insert((*m_perfectLd).end(), perfectLd.begin(), perfectLd.end());
    LinkageThread::mtx.unlock();

}


void *LinkageThread::triangularProcess(void *in){
    ((LinkageThread *) in)->triangularProcess();
    return nullptr;
}

void *LinkageThread::rectangularProcess(void *in){
    ((LinkageThread *) in)->rectangularProcess();
    return nullptr;
}


