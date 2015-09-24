#include "linkage.h"

Linkage::Linkage(size_t thread):m_thread(thread){}


void Linkage::Initialize(boost::ptr_deque<Genotype> &genotype, const size_t &prevResiduals, const boost::ptr_vector<Interval> &blockSize){
	if(genotype.empty()){
        std::runtime_error("Cannot build LD without genotypes")
	}
	//Use previous informations
    if(prevResiduals == 0){
        m_linkage = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
        m_linkageSqrt = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
    }
    else{
        blockSize
		Eigen::MatrixXd temp = m_linkage.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
        //Currently we only need the R2 matrix, so we will ignore the R matrix
		//temp = m_linkageSqrt.bottomRightCorner(prevResiduals,prevResiduals);
		//m_linkageSqrt= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		//m_linkageSqrt.topLeftCorner(prevResiduals, prevResiduals) = temp;

    }
}


void Linkage::Construct(boost::ptr_deque<Genotype> &genotype, const size_t &genotypeIndex, const size_t&remainedLD, const boost::ptr_vector<Interval> &blockSize, bool correction){
	if(genotype.empty())    throw std::runtime_error("There is no genotype to work on");
    for(size_t i = 0; i < m_thread; ++i){

    }
}
