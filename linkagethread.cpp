#include "linkagethread.h"
std::mutex LinkageThread::mtx;

void LinkageThread::Addstart(size_t i){
	m_startLoc.push_back(i);
}

LinkageThread::LinkageThread(bool correction, const size_t blockEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList):m_correction(correction), m_boundEnd(blockEnd), m_ldMatrix(ldMatrix), m_genotype(genotype), m_snpLoc(snpLoc), m_perfectLd(perfectLd), m_snpList(snpList){}

LinkageThread::LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList):m_correction(correction), m_snpStart(snpStart), m_snpEnd(snpEnd), m_boundStart(boundStart), m_boundEnd(boundEnd), m_ldMatrix(ldMatrix), m_genotype(genotype), m_snpLoc(snpLoc), m_perfectLd(perfectLd), m_snpList(snpList){}

LinkageThread::~LinkageThread(){}

void LinkageThread::triangularProcess(){
    size_t sizeOfItem = m_startLoc.size();
    std::vector<size_t> perfectLd;
    for(size_t i = 0; i < sizeOfItem; ++i){
        size_t j = m_startLoc[i];
        (*m_ldMatrix)(j,j) = 1.0;
        for(size_t k = m_boundEnd-1; k > j; --k){
            if((*m_ldMatrix)(j,k) == 0.0){
                double rSquare=(*m_genotype)[j]->Getr( (*m_genotype)[k], m_correction);
                if(std::fabs(rSquare-1.0) < G_EPSILON_DBL ){
                    perfectLd.push_back(k);
                    LinkageThread::mtx.lock();
                        (*m_snpList)[(*m_snpLoc)[k]]->shareHeritability((*m_snpList)[(*m_snpLoc)[j]]);
                    LinkageThread::mtx.unlock();
                }
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
            }
            else if((*m_ldMatrix)(i,j) == 0.0){
                double rSquare = (*m_genotype)[i]->Getr((*m_genotype)[j],m_correction);
                if(std::fabs(rSquare-1.0) < G_EPSILON_DBL){
                    perfectLd.push_back(j);
                    LinkageThread::mtx.lock();
                        (*m_snpList)[(*m_snpLoc)[j]]->shareHeritability((*m_snpList)[(*m_snpLoc)[i]]);
                    LinkageThread::mtx.unlock();
                }
				(*m_ldMatrix)(i,j) = rSquare;
				(*m_ldMatrix)(j,i) = rSquare;

            }
            else break;
        }
    }
    LinkageThread::mtx.lock();//DEBUG
    (*m_perfectLd).insert((*m_perfectLd).end(), perfectLd.begin(), perfectLd.end());
    LinkageThread::mtx.unlock();

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


