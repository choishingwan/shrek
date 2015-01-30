#include "linkagethread.h"

std::mutex LinkageThread::mtx;

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

void LinkageThread::Seteffective(double effective){
    *m_effectiveNumber += effective;
}

LinkageThread::LinkageThread(bool correction, const size_t blockEnd, double *effectiveNum, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype):m_correction(correction), m_boundEnd(blockEnd), m_effectiveNumber(effectiveNum), m_ldMatrix(ldMatrix), m_genotype(genotype){}

LinkageThread::LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, double *effectiveNum, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype):m_correction(correction), m_snpStart(snpStart), m_snpEnd(snpEnd), m_boundStart(boundStart), m_boundEnd(boundEnd), m_effectiveNumber(effectiveNum), m_ldMatrix(ldMatrix), m_genotype(genotype){}

LinkageThread::~LinkageThread()
{
	//dtor
}


void *LinkageThread::triangularProcess(void *in){
    mtx.lock();
    struct LinkageThread *input = (LinkageThread *) in;
    bool correction = input->Getcorrection();
    size_t endOfProcess = input->GetboundEnd();
    Eigen::MatrixXd *matrix = input->Getld();
    std::deque<Genotype*> *geno = input->Getgenotype();

    size_t sizeOfItem = input->GetsizeOfStart();
    double effective=0.0;
    for(size_t i = 0; i < sizeOfItem; ++i){
        size_t j = input->GetstartLoc(i);
        (*matrix)(j,j) = 1.0;
        effective+= 1.0;
        for(size_t k = j+1; k < endOfProcess; ++k){
            double rSquare=(*geno)[j]->Getr( (*geno)[k], correction);
            (*matrix)(j, k) = rSquare;
            (*matrix)(k,j) = rSquare;
            effective += 2.0*rSquare;
        }
    }
    mtx.unlock();

    mtx.lock();
	input->Seteffective(effective);
    mtx.unlock();
    return nullptr;
}



void *LinkageThread::rectangularProcess(void *in){
    struct LinkageThread *input = (LinkageThread *) in;

    bool correction = input->Getcorrection();
    size_t snpStart = input->GetsnpStart();
    size_t snpEnd = input->GetsnpEnd();
    size_t boundStart = input->GetboundStart();
    size_t boundEnd = input->GetboundEnd();
    Eigen::MatrixXd *matrix = input->Getld();
    std::deque<Genotype* > *geno = input->Getgenotype();
    double effective = 0.0;

    for(size_t i = snpStart; i < snpEnd; ++i){
        for(size_t j = boundStart; j < boundEnd; ++j){
            if(j == i) (*matrix)(i,i) = 1.0; //Let's just assume that it is duplicated
            else{
                double rSquare = (*geno)[i]->Getr((*geno)[j],correction);
				(*matrix)(i,j) = rSquare;
				(*matrix)(j,i) = rSquare;
				effective += 2.0*rSquare;
            }
        }
    }

	mtx.lock();
    input->Seteffective(effective);
    mtx.unlock();
    return nullptr;
}


