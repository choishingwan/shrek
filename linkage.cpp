#include "linkage.h"

Linkage::Linkage(size_t thread):m_thread(thread){}
Linkage::~Linkage(){}
std::mutex Linkage::mtx;


void Linkage::Initialize(const boost::ptr_deque<Genotype> &genotype, const size_t &prevResiduals){
	if(genotype.empty()){
        std::runtime_error("Cannot build LD without genotypes");
	}
	//Use previous informations
    if(prevResiduals == 0){
        m_linkage = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
        //m_linkageSqrt = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
    }
    else{
		Eigen::MatrixXd temp = m_linkage.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
        //Currently we only need the R2 matrix, so we will ignore the R matrix
		//temp = m_linkageSqrt.bottomRightCorner(prevResiduals,prevResiduals);
		//m_linkageSqrt= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		//m_linkageSqrt.topLeftCorner(prevResiduals, prevResiduals) = temp;

    }
//    std::ofstream debug("DEBUG");
//    debug << m_linkage << std::endl;
//    debug.close();
}

void Linkage::buildLd(const bool correction, const size_t vStart, const size_t vEnd, const size_t hEnd, const boost::ptr_deque<Genotype> &genotype, const std::deque<size_t> &ldLoc){
    //Will work on all SNPs from vStart to hEnd
    size_t genotypeStart = genotype.size();
    size_t genotypeEnd = genotype.size(); //Bound
    size_t bottom = genotype.size(); //Bound
    for(size_t i = 0; i < ldLoc.size(); ++i){
        if(ldLoc[i]==vStart) genotypeStart = i;
        else if(ldLoc[i] == vEnd) bottom = i+1;
        if(ldLoc[i]==hEnd){
            genotypeEnd = i+1;
            break;
        }
    }

//    Linkage::mtx.lock();
//    std::cerr << "Block info: " << genotypeStart << "\t" << genotypeEnd << "\t" << bottom << std::endl;
//    Linkage::mtx.unlock();
//
    for(size_t i = genotypeStart; i < bottom; ++i){
        m_linkage(i,i) = 1.0;
        for(size_t j = i+1; j< genotypeEnd; ++j){
            if(m_linkage(i,j) == 0.0){ //Again, this is problematic because you are finding equal of double
                double rSquare=0.0;
                double r =0.0;
                genotype[i].GetbothR(&(genotype[j]),correction, r, rSquare);
                m_linkage(i,j) = rSquare;
                m_linkage(j,i) = rSquare;
            }
        }
    }
}

void Linkage::Construct(const boost::ptr_deque<Genotype> &genotype, const size_t &genotypeIndex, const size_t& remainedLD, const boost::ptr_vector<Interval> &blockInfo, const bool correction, const std::deque<size_t> &ldLoc){
	if(genotype.empty())    throw std::runtime_error("There is no genotype to work on");
    size_t startRange =  genotypeIndex;
    size_t endRange=0;
    size_t range = m_thread;
    if(remainedLD==0)range +=2;
    else{
        //Something was left behind, and they must be 2 blocks before
        startRange -= 2;//Two steps is only right if I have remove stuffs
    }
    std::string currentChr = blockInfo[startRange].getChr();
    size_t boundHunter = genotypeIndex;
    for(; boundHunter< genotypeIndex+range && boundHunter < blockInfo.size(); ++boundHunter){
            if(blockInfo[boundHunter].getChr().compare(currentChr)!=0){
            boundHunter= boundHunter-1; //This mean we working on the last block of this chromosome
            break;
        }
    }
    endRange = boundHunter;
    // So now the startRange and endRange will contain the index of the intervals to include in the LD construction
    // Each thread will process one sausage       \------------|
    //                                             \-----------|

    //Change this into threading
    std::vector<std::thread> threadStore;
    //Launch a group of threads
    size_t threadRunCounter = startRange;
    while(threadRunCounter < endRange){
        while(threadStore.size() < m_thread && threadRunCounter < endRange){ //On purposely leave 1 thread out for the main
            if(threadRunCounter+2 < endRange){
                threadStore.push_back(std::thread(&Linkage::buildLd, this, correction, blockInfo[threadRunCounter].getStart(), blockInfo[threadRunCounter].getEnd()-1,blockInfo[threadRunCounter+2].getEnd(), std::cref(genotype), std::cref(ldLoc)));
            }
            else if(threadRunCounter+1 < endRange){
                threadStore.push_back(std::thread(&Linkage::buildLd, this, correction, blockInfo[threadRunCounter].getStart(), blockInfo[threadRunCounter].getEnd()-1,blockInfo[threadRunCounter+1].getEnd(), std::cref(genotype), std::cref(ldLoc)));
            }
            else{
                threadStore.push_back(std::thread(&Linkage::buildLd, this, correction, blockInfo[threadRunCounter].getStart(), blockInfo[threadRunCounter].getEnd(),blockInfo[threadRunCounter].getEnd(), std::cref(genotype), std::cref(ldLoc)));
            }
            threadRunCounter++;
        }

        for (size_t j = 0; j < threadStore.size(); ++j) {
            threadStore[j].join();
        }
        threadStore.clear();
    }
}

void Linkage::print(){
    std::cout << m_linkage << std::endl;
}

void Linkage::solve(const size_t loc, const size_t length, const Eigen::MatrixXd &betaEstimate, Eigen::MatrixXd &heritability, Eigen::MatrixXd &effectiveNumber, Eigen::VectorXd &ldScore) const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m_linkage.block(loc, loc, length, length));
    double tolerance = std::numeric_limits<double>::epsilon() * length * es.eigenvalues().array().maxCoeff();
    Eigen::MatrixXd rInverse = es.eigenvectors()*(es.eigenvalues().array() > tolerance).select(es.eigenvalues().array().inverse(), 0).matrix().asDiagonal() * es.eigenvectors().transpose();
    /** Calculate the h vector here **/
    heritability= rInverse*betaEstimate.block(loc, 0,length,betaEstimate.cols());
    Eigen::MatrixXd error =m_linkage.block(loc, loc, length, length)*heritability - betaEstimate.block(loc, 0,length,betaEstimate.cols());
	double bNorm = betaEstimate.block(loc, 0,length,betaEstimate.cols()).norm();
    double relative_error = error.norm() / bNorm;
    double prev_error = relative_error+1;
    Eigen::MatrixXd update = heritability;
    //Iterative improvement, arbitrary max iteration
    size_t maxIter = 300;
    size_t iterCount = 0;
    while(relative_error < prev_error && iterCount < maxIter){
        prev_error = relative_error;
        update.noalias()=rInverse*(-error);
        relative_error = 0.0;
        error.noalias()= m_linkage.block(loc, loc, length, length)*(heritability+update) - betaEstimate.block(loc, 0,length,betaEstimate.cols());
        relative_error = error.norm() / bNorm;
        if(relative_error < 1e-300) relative_error = 0;
        heritability = heritability+update;
        iterCount++;
    }


    Eigen::MatrixXd vecOfOne = Eigen::MatrixXd::Constant(length,betaEstimate.cols(), 1.0);
    double eNorm = vecOfOne.norm();
    /** Calculate the effective number here **/
    effectiveNumber = rInverse*vecOfOne;
    error = m_linkage.block(loc, loc, length, length)*effectiveNumber-vecOfOne;
    relative_error = error.norm()/eNorm;
    prev_error = relative_error+1;
    update=effectiveNumber;
    iterCount = 0;
    while(relative_error < prev_error && iterCount < maxIter){
        prev_error = relative_error;
        update.noalias()=rInverse*(-error);
        relative_error = 0.0;
        error.noalias()= m_linkage.block(loc, loc, length, length)*(effectiveNumber+update) - vecOfOne;
        relative_error = error.norm() / eNorm;
        if(relative_error < 1e-300) relative_error = 0;
        effectiveNumber = effectiveNumber+update;
        iterCount++;
    }
    /** Calculate LDSCore **/
    ldScore =m_linkage.block(loc, loc, length, length).colwise().sum();
}
