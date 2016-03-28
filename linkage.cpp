#include "linkage.h"

std::mutex Linkage::linkageMtx;
size_t Linkage::check = 0;
Linkage::Linkage(size_t thread):m_thread(thread){}
Linkage::~Linkage(){}

void Linkage::print(){ std::cout << m_linkage << std::endl; }
void Linkage::print(size_t start, size_t ending, std::string name){
    std::ofstream checking;
    checking.open(name.c_str());
    checking << m_linkage.submat(start,start,ending,ending) << std::endl;
    checking.close();
    std::string inv = name+".rinv";
    checking.open(inv.c_str());
    checking << m_rInv << std::endl;
    checking.close();
}

void Linkage::computeLd(const boost::ptr_deque<Genotype> &genotype, const std::deque<size_t> &snpLoc, size_t verticalStartIndex, size_t verticalEndBound, size_t horizontalStartIndex, boost::ptr_vector<Snp> &snpList, const bool &correction, std::vector<size_t> &perfectLd){
    std::vector<size_t> perfectLDBufferStore; //This stores the SNP that are going to be removed because of perfect LD
    std::vector<size_t> perfectLDRemainStore; //This stores the survivor from the perfect LD removal
    size_t horizontalEndBound = snpLoc.size();
    for(size_t i = verticalStartIndex; i < verticalEndBound; ++i){
        size_t horizontalStart = (i>horizontalStartIndex)? i:horizontalStartIndex;
        for(size_t j = horizontalStart; j < horizontalEndBound; ++j){
                if(m_linkage(i,j)==0.0){  // only calculate the R2 when this is a new spot
                    if(i==j){
                        m_linkage(i,j) = 1.0;
                        m_linkageSqrt(i,j) = 1.0;
                    }
                    else{
                        double r = 0.0;
                        double r2 = 0.0;
                        genotype.at(i).GetbothR(genotype.at(j), correction, r, r2);
                        m_linkage(i,j) = r2;
                        m_linkageSqrt(i,j) = r;
                        m_linkage(j,i) = r2;
                        m_linkageSqrt(j,i) = r;
                        if(std::fabs(r-1.0) < 1e-6 || r > 1.0 || r < -1.0){ //Perfect LD if R = 1 or > abs(1)
                            perfectLDBufferStore.push_back(j);
                            perfectLDRemainStore.push_back(i);
                        }
                    }
                }
        }
    }
//  Now perform the thread safe recording of the perfect LD and update the beta
    linkageMtx.lock();
        for(size_t i = 0; i < perfectLDBufferStore.size(); ++i) snpList[snpLoc[perfectLDRemainStore[i]]].shareHeritability(snpList[snpLoc[perfectLDBufferStore[i]]]);
        perfectLDRemainStore.clear();
        perfectLd.insert(perfectLd.end(), perfectLDBufferStore.begin(), perfectLDBufferStore.end());
    linkageMtx.unlock();

}

void Linkage::construct(boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::vector<size_t> &boundary, boost::ptr_vector<Snp> &snpList, const bool correction,bool &boundCheck){
    // Here we try to prepare for the calculation of LD
    // This function is mainly responsible for the thread distribution
    if(genotype.empty()) throw std::runtime_error("Cannot build LD without genotypes");
    // If this is the first of its kind, then we need to initialize the matrix
    size_t genoSize = genotype.size();
	if(m_linkage.n_cols==0){
        m_linkage = arma::mat(genoSize, genoSize,arma::fill::eye);
        m_linkageSqrt = arma::mat(genoSize, genoSize, arma::fill::eye);
    }
    else{ // Otherwise, just make sure the matrix is of the right size
        m_linkage.resize(genoSize, genoSize);
        m_linkageSqrt.resize(genoSize, genoSize);
    }
    std::vector<size_t> perfectLd; // This is the vector use to store the perfect LD information
    int boundSize = boundary.size();
    size_t startBoundIndex = (boundSize>3)? boundary[boundSize-3]: boundary.front();
    size_t snpLocSize = snpLoc.size();
    int numSnpRequired = snpLocSize - startBoundIndex;
    assert(numSnpRequired>0 && "Must have non-zero number of SNPs to work on");
    std::vector<std::thread> threadStore;
    int jobBlockSize  = numSnpRequired /m_thread;
    int remainSnps = numSnpRequired%m_thread;
    int currentBlock = startBoundIndex;
    for(size_t i = 0; i < m_thread && currentBlock < snpLocSize; ++i){
        if(remainSnps>0){
            threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), currentBlock, currentBlock+jobBlockSize+1, boundary.back(), std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
            --remainSnps;
            currentBlock+=jobBlockSize+1;
        }
        else{
            threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), currentBlock, currentBlock+jobBlockSize, boundary.back(), std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
            currentBlock+=jobBlockSize;
        }
    }
    for (size_t i = 0; i < threadStore.size(); ++i) threadStore[i].join();
    threadStore.clear();
    std::sort(perfectLd.begin(), perfectLd.end());
    perfectLd.erase( std::unique( perfectLd.begin(), perfectLd.end() ), perfectLd.end() );
    if(perfectLd.size()!=0) perfectRemove(perfectLd, genotype, snpLoc, boundary, snpList, boundCheck);
}


void Linkage::perfectRemove(std::vector<size_t> &perfectLd, boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::vector<size_t> &boundary, boost::ptr_vector<Snp> &snpList, bool &boundCheck){
    // First, update the matrix
    size_t genoSize = genotype.size();
    size_t cI = 0;
    size_t perfectIndexI = 0;
    //This script can actually be performed using multi-thread. But we will leave this out for now
    for(size_t i = 0; i < genoSize; ++i){
        if(i == perfectLd[perfectIndexI]) perfectIndexI++;// This is perfect LD, we will skip it
        else{ // The current I is required, go over the remaining
            size_t cJ = cI, perfectIndexJ = perfectIndexI; // Because of symmetric matrix, we only start adding at cI
            for(size_t j = i; j < genoSize; ++j){
                if(j == perfectLd[perfectIndexJ]) perfectIndexJ++;
                else{ //This is required
                    m_linkage(cJ,cI) = m_linkage(i,j);
                    m_linkageSqrt(cJ,cI)= m_linkageSqrt(i,j);;
                    m_linkage(cI,cJ) = m_linkage(i,j);
                    m_linkageSqrt(cI,cJ) = m_linkageSqrt(i,j);
                    cJ++;
                }
            }
            cI++;
        }
    }
    size_t numPerfect = perfectLd.size();
    // Now resize the matrix to remove unwanted stuff
//    m_linkage.resize(genoSize-numPerfect, genoSize-numPerfect);
//    m_linkageSqrt.resize(genoSize-numPerfect, genoSize-numPerfect);
    m_linkage = m_linkage.submat(0,0,cI,cI); //Feel like this is safer
    m_linkageSqrt =m_linkageSqrt.submat(0,0,cI,cI);
    // Now update the snpLoc and genotype
    for(size_t i = 0; i < numPerfect; ++i){
        snpLoc.erase(snpLoc.begin()+(perfectLd[i]-i));
        genotype.erase(genotype.begin()+(perfectLd[i]-i));
    }
    // Now carefully check if the last boundary is in perfect LD
    if(boundary.size() > 1){
        size_t snpLocSize = snpLoc.size();
        size_t smaller = 0;
        for(size_t i = 0; i < perfectLd.size(); ++i){
            if(perfectLd[i]==boundary.back()){
                boundCheck=true;
                boundary.back()++;
            }
            if(perfectLd[i] < boundary.back()) smaller++;
        }
        // Now change the boundary
        boundary.back()-= smaller;
    }
}

void Linkage::clear(){
    m_linkage.clear();
    m_linkageSqrt.clear();
}

void Linkage::clear(size_t nRemoveElements){
    size_t lastIndex = m_linkage.n_cols-1;
    m_linkage = m_linkage.submat(nRemoveElements, nRemoveElements,lastIndex,lastIndex);
    m_linkageSqrt = m_linkageSqrt.submat(nRemoveElements, nRemoveElements,lastIndex,lastIndex);
}

void Linkage::complexSE(size_t startIndex, size_t endIndex, const arma::vec &nSample, const arma::vec &tStat, arma::mat &varResult){
    varResult = m_rInv*arma::diagmat(nSample)*(4*m_linkageSqrt.submat(startIndex,startIndex,endIndex,endIndex)%(tStat*tStat.t())-2*m_linkage.submat(startIndex,startIndex,endIndex, endIndex))*arma::diagmat(nSample)*m_rInv;
    m_rInv.clear();
}

void Linkage::effectiveSE(size_t startIndex, size_t endIndex, arma::vec &varResult){
    if(startIndex > m_linkage.n_cols) throw "Start coordinates exceeds the matrix size";
    arma::vec effectiveNumber(endIndex-startIndex, arma::fill::ones);
    varResult = m_rInv*effectiveNumber;
    arma::vec errorVec = m_linkage.submat( startIndex, startIndex, endIndex, endIndex )*varResult - effectiveNumber;
 	double oriNorm = norm(effectiveNumber);
    double relative_error = norm(errorVec) / oriNorm;
    double prev_error = relative_error+1;
    arma::vec updateVec;
    int iterCount = 0, maxIter = 10000;
    while(relative_error < prev_error && iterCount < maxIter){
        prev_error = relative_error;
        updateVec=m_rInv*(-errorVec);
        errorVec= m_linkage.submat( startIndex, startIndex, endIndex, endIndex )*(varResult+updateVec) - effectiveNumber;
        relative_error = norm(errorVec) / oriNorm;
        if(relative_error < 1e-300) relative_error = 0; // 1e-300 is more than enough...
        varResult = varResult+updateVec;
        iterCount++; // This is to avoid infinite looping
    }
    m_rInv.clear();
}


//double Linkage::m_tolerance = 0.0;
void Linkage::decompose(size_t start, const arma::vec &fStat, arma::vec &heritResult){
    if(start > m_linkage.n_cols) throw "Start coordinates exceeds the matrix size";
    size_t endOfBlock = start+fStat.n_elem-1;
    m_rInv=arma::pinv((arma::mat)m_linkage.submat( start, start, endOfBlock, endOfBlock ));
//    arma::vec eigval, eigvec;
//    arma::eig_sym( eigval, eigvec, (arma::mat)m_linkage.submat( start, start, endOfBlock, endOfBlock ) );
//    arma::vec eigval;
//    arma::mat U, V;
//    arma::svd(U,eigval,V,subMat);
//    Linkage::m_tolerance = std::numeric_limits<double>::epsilon() * fStat.n_elem *arma::max(eigval);
//    eigval.for_each( [](arma::mat::elem_type& val) { if(val < Linkage::m_tolerance){val=0;}else{val = 1/val; }} );
//    m_rInv = eigvec*arma::diagmat(eigval) * eigvec.t();
//    m_rInv = V*arma::diagmat(eigval) * U.t();
    heritResult = m_rInv * fStat;
    /** Calculate the h vector here **/
    arma::mat error = m_linkage.submat( start, start, endOfBlock, endOfBlock )*heritResult - fStat;
 	double oriNorm = norm(fStat);
    double relative_error = norm(error) / oriNorm;
    double prev_error = relative_error+1;
    arma::mat update;
    int iterCount = 0, maxIter = 10000;
    while(relative_error < prev_error && iterCount < maxIter){
        prev_error = relative_error;
        update=m_rInv*(-error);
        error= m_linkage.submat( start, start, endOfBlock, endOfBlock )*(heritResult+update) - fStat;
        relative_error = norm(error) / oriNorm;
        if(relative_error < 1e-300) relative_error = 0; // 1e-300 is more than enough...
        heritResult = heritResult+update;
        iterCount++; // This is to avoid infinite looping
    }

}

