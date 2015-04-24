#include "linkagethread.h"
std::mutex LinkageThread::mtx;

void LinkageThread::Addstart(size_t i){ m_startLoc.push_back(i); }

LinkageThread::LinkageThread(bool correction, const size_t blockEnd, Eigen::MatrixXd *ldMatrix, Eigen::MatrixXd *varMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList):m_correction(correction), m_boundEnd(blockEnd), m_ldMatrix(ldMatrix), m_varLdMatrix(varMatrix),m_genotype(genotype), m_snpLoc(snpLoc), m_perfectLd(perfectLd), m_snpList(snpList){}

LinkageThread::LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, Eigen::MatrixXd *ldMatrix, Eigen::MatrixXd *varMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList):m_correction(correction), m_snpStart(snpStart), m_snpEnd(snpEnd), m_boundStart(boundStart), m_boundEnd(boundEnd), m_ldMatrix(ldMatrix), m_varLdMatrix(varMatrix), m_genotype(genotype), m_snpLoc(snpLoc), m_perfectLd(perfectLd), m_snpList(snpList){}

LinkageThread::~LinkageThread(){}

void LinkageThread::triangularProcess(){
    size_t sizeOfItem = m_startLoc.size();
    std::vector<size_t> perfectLd;
    for(size_t i = 0; i < sizeOfItem; ++i){
        size_t j = m_startLoc[i];
        (*m_ldMatrix)(j,j) = 1.0;
        (*m_varLdMatrix)(j,j) = Linkage::VarianceR2(1.0,(*m_genotype)[j]->GetnumSample(),0);
        for(size_t k = m_boundEnd-1; k > j; --k){
            if((*m_ldMatrix)(j,k) == 0.0){
                size_t numSample=0;
                double rSquare=(*m_genotype)[j]->GetrSq( (*m_genotype)[k], m_correction, numSample);
                if(std::fabs(rSquare-1.0) < G_EPSILON_DBL ){
                    perfectLd.push_back(k);
                    LinkageThread::mtx.lock();
                        (*m_snpList)[(*m_snpLoc)[k]]->shareHeritability((*m_snpList)[(*m_snpLoc)[j]]);
                    LinkageThread::mtx.unlock();
                }
                (*m_varLdMatrix)(j,k) =Linkage::VarianceR2(rSquare, numSample,1);
                (*m_varLdMatrix)(k,j) = (*m_varLdMatrix)(j,k);
                (*m_ldMatrix)(j, k) = rSquare;
                (*m_ldMatrix)(k,j) = rSquare;
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
                (*m_varLdMatrix)(i,i) = Linkage::VarianceR2(1.0,(*m_genotype)[i]->GetnumSample(),0);
            }
            else if((*m_ldMatrix)(i,j) == 0.0){
                size_t numSample = 0;
                double rSquare = (*m_genotype)[i]->GetrSq((*m_genotype)[j],m_correction, numSample);
                if(std::fabs(rSquare-1.0) < G_EPSILON_DBL){
                    perfectLd.push_back(j);
                    LinkageThread::mtx.lock();
                        (*m_snpList)[(*m_snpLoc)[j]]->shareHeritability((*m_snpList)[(*m_snpLoc)[i]]);
                    LinkageThread::mtx.unlock();
                }
                (*m_varLdMatrix)(j,i) = Linkage::VarianceR2(rSquare, numSample,1);
                (*m_varLdMatrix)(i,j) = (*m_varLdMatrix)(j,i);
				(*m_ldMatrix)(i,j) = rSquare;
				(*m_ldMatrix)(j,i) = rSquare;

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


