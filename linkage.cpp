#include "linkage.h"

std::mutex Linkage::linkageMtx;

Linkage::Linkage(size_t thread, size_t blockSize):m_thread(thread), m_blockSize(blockSize){}
Linkage::~Linkage(){}

void Linkage::print(){ std::cout << m_linkage << std::endl; }

//void Linkage::computeLd(const boost::ptr_list<Genotype> &genotype, const std::list<size_t> &snpLoc, size_t startIndex, size_t verEnd, size_t horistart, size_t horiEnd, boost::ptr_vector<Snp> &snpList, const bool &correction, std::vector<size_t> &perfectLd){
void Linkage::computeLd(const boost::ptr_deque<Genotype> &genotype, const std::deque<size_t> &snpLoc, size_t startIndex, size_t verEnd, size_t horistart, boost::ptr_vector<Snp> &snpList, const bool &correction, std::vector<size_t> &perfectLd){
    std::vector<size_t> perfectBuff;
    size_t horiEnd = snpLoc.size();


    for(size_t i = startIndex; i < verEnd; ++i){
        size_t nextStart = (i>horistart)? i:horistart;
        size_t j = nextStart;
        for(size_t j = nextStart; j < horiEnd; ++j){
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
                        // Now check if it is perfect LD
                        // Store the index j into perfectBuff
                        if(std::fabs(r-1.0) < 1e-6 || r > 1.0 || r < -1.0){
                            perfectBuff.push_back(j);
                            snpList.at(snpLoc.at(i)).shareHeritability(snpList.at(snpLoc.at(j)));
                        }
                    }
                }
        }
    }
    // Now we can insert the buff back to the result
    linkageMtx.lock();
        perfectLd.insert(perfectLd.end(), perfectBuff.begin(), perfectBuff.end());
    linkageMtx.unlock();

}

void Linkage::construct(boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::deque<size_t> &boundary, boost::ptr_vector<Snp> &snpList, const bool correction,bool &boundCheck){
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
    size_t startBoundIndex = (boundSize-3>0)? boundary[boundSize-3]: boundary.front();
    size_t snpLocSize = snpLoc.size();
    int numSnpRequired = snpLocSize - startBoundIndex;
    assert(numSnpRequired>0 && "Must have non-zero number of SNPs to work on");
    std::vector<std::thread> threadStore;
    if(numSnpRequired < m_thread){
        // We will use as much threads as possible
        for(size_t i = startBoundIndex; i < snpLocSize; ++i){
            threadStore.push_back(
                std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc),
                i, i+1, boundary.back(),  // Vertical start, vertical end, horizontal start
                std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
        }
    }
    else{
        // We actually need to divide the job
        int jobBlock = numSnpRequired/m_thread;
        int remainder = numSnpRequired%m_thread;
        int currentBlock = startBoundIndex;
        for(size_t i = 0; i < m_thread; ++i){
            if(remainder>0){
                threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), currentBlock, currentBlock+jobBlock+1, boundary.back(), std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
                --remainder;
                currentBlock+=jobBlock+1;
            }
            else{
                threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), currentBlock, currentBlock+jobBlock, boundary.back(), std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
                currentBlock+=jobBlock;
            }
        }
    }
    size_t threadSize = threadStore.size();
    for (size_t i = 0; i < threadSize; ++i) threadStore[i].join();
    threadStore.clear();

    // We sort the index for easier management (they are not sorted because of multi threading
    std::sort(perfectLd.begin(), perfectLd.end());
    // Remove duplicated index
    perfectLd.erase( std::unique( perfectLd.begin(), perfectLd.end() ), perfectLd.end() );

    // If something are in perfect LD, perform the perfect LD removal
    if(perfectLd.size()!=0) perfectRemove(perfectLd, genotype, snpLoc, boundary, snpList, boundCheck);

}



void Linkage::perfectRemove(std::vector<size_t> &perfectLd, boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::deque<size_t> &boundary, boost::ptr_vector<Snp> &snpList, bool &boundCheck){
    // First, update the matrix
    size_t genoSize = genotype.size();
    size_t cI = 0;
    size_t perfectIndexI = 0;
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
    m_linkage.resize(genoSize-numPerfect, genoSize-numPerfect);
    m_linkageSqrt.resize(genoSize-numPerfect, genoSize-numPerfect);
    // Now update the snpLoc and genotype
    size_t itemRemoved = 0;
    for(size_t i = 0; i < numPerfect; ++i){
        snpLoc.erase(snpLoc.begin()+(perfectLd[i]-itemRemoved));
        genotype.erase(genotype.begin()+(perfectLd[i]-itemRemoved));
        itemRemoved++;
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

    size_t blockSize = m_linkage.n_cols-1;
    m_linkage = m_linkage.submat(nRemoveElements, nRemoveElements,blockSize,blockSize);
    m_linkageSqrt = m_linkageSqrt.submat(nRemoveElements, nRemoveElements,blockSize,blockSize);
}

void Linkage::decompose(size_t start, const arma::vec &fStat, arma::vec &heritResult, arma::vec &varResult){
    if(start > m_linkage.n_cols) throw "Start coordinates exceeds the matrix size";
    size_t endOfBlock = fStat.n_elem-1;
    arma::mat rInv=pinv((arma::mat)m_linkage.submat( start, start, endOfBlock,endOfBlock ));
    heritResult = rInv * fStat;
    arma::mat error = m_linkage.submat( start, start, endOfBlock, endOfBlock )*heritResult - fStat;
 	double oriNorm = norm(fStat);
    double relative_error = norm(error) / oriNorm;
    double prev_error = relative_error+1;
    arma::mat update;
    int iterCount = 0, maxIter = 100;
    while(relative_error < prev_error && iterCount < maxIter){
        prev_error = relative_error;
        update=rInv*(-error);
        error= m_linkage.submat( start, start, endOfBlock, endOfBlock )*(heritResult+update) - fStat;
        relative_error = norm(error) / oriNorm;
        if(relative_error < 1e-300) relative_error = 0; // 1e-300 is more than enough...
        heritResult = heritResult+update;
        iterCount++; // This is to avoid infinite looping
    }
    // Now calculate the variance
    arma::vec effectiveNumber(heritResult.n_elem, arma::fill::ones);
    varResult = rInv*effectiveNumber;
    arma::vec errorVec = m_linkage.submat( start, start, endOfBlock, endOfBlock )*varResult - effectiveNumber;
 	oriNorm = norm(effectiveNumber);
    relative_error = norm(errorVec) / oriNorm;
    prev_error = relative_error+1;
    arma::vec updateVec;
    iterCount = 0;
    while(relative_error < prev_error && iterCount < maxIter){
        prev_error = relative_error;
        updateVec=rInv*(-errorVec);
        errorVec= m_linkage.submat( start, start, endOfBlock, endOfBlock )*(varResult+updateVec) - effectiveNumber;
        relative_error = norm(errorVec) / oriNorm;
        if(relative_error < 1e-300) relative_error = 0; // 1e-300 is more than enough...
        varResult = varResult+updateVec;
        iterCount++; // This is to avoid infinite looping
    }

}

void Linkage::decompose(size_t start, const arma::vec &zStat, const arma::vec &fStat, const arma::vec &nSample, arma::vec &heritResult, arma::mat &varResult){
    if(start > m_linkage.n_cols) throw "Start coordinates exceeds the matrix size";
    size_t endOfBlock = start+fStat.n_elem-1;
    arma::mat rInv=pinv((arma::mat)m_linkage.submat( start, start, endOfBlock, endOfBlock ));
    heritResult = rInv * fStat;
    arma::mat error = m_linkage.submat( start, start, endOfBlock, endOfBlock )*heritResult - fStat;
 	double oriNorm = norm(fStat);
    double relative_error = norm(error) / oriNorm;
    double prev_error = relative_error+1;
    arma::mat update;
    int iterCount = 0, maxIter = 100;
    while(relative_error < prev_error && iterCount < maxIter){
        prev_error = relative_error;
        update=rInv*(-error);
        error= m_linkage.submat( start, start, endOfBlock, endOfBlock )*(heritResult+update) - fStat;
        relative_error = norm(error) / oriNorm;
        if(relative_error < 1e-300) relative_error = 0; // 1e-300 is more than enough...
        heritResult = heritResult+update;
        iterCount++; // This is to avoid infinite looping
    }

    // Now calculate the varaince, this is for the complicated variance
    varResult = rInv*arma::diagmat(nSample)*(4*m_linkageSqrt.submat(start,start,endOfBlock,endOfBlock)%(zStat*zStat.t())-2*m_linkage.submat(start,start,endOfBlock, endOfBlock))*arma::diagmat(nSample)*rInv;


    //This part is here to safe guard any stupid mistakes I made
    if(isnan(arma::accu(varResult))){
        std::cout << varResult << std::endl;
        std::cout << "sample:"  << std::endl;
        std::cout << arma::diagmat(nSample) << std::endl;
        std::cout << "R2" << std::endl;
        std::cout << m_linkage.submat(start,start,endOfBlock, endOfBlock) << std::endl;
        std::cout << "R" << std::endl;
        std::cout << m_linkageSqrt.submat(start,start,endOfBlock,endOfBlock) << std::endl;
        std::cout << "Z" << std::endl;
        std::cout << (zStat*zStat.t()) << std::endl;
        std::cout << "Inverse" << std::endl;
        std::cout << rInv << std::endl;
        std::cout << "Fstat"<< std::endl;
        std::cout << fStat << std::endl;
    }
}

void Linkage::decompose(size_t start, const arma::vec &zStat, const arma::vec &fStat, const arma::vec &nSample, arma::vec &heritResult, arma::mat &varResult, arma::mat &addVarResult){
    if(start > m_linkage.n_cols) throw "Start coordinates exceeds the matrix size";
    size_t endOfBlock = start+fStat.n_elem-1;
    arma::mat rInv=pinv((arma::mat)m_linkage.submat( start, start, endOfBlock, endOfBlock ));
    heritResult = rInv * fStat;
    arma::mat error = m_linkage.submat( start, start, endOfBlock, endOfBlock )*heritResult - fStat;
 	double oriNorm = norm(fStat);
    double relative_error = norm(error) / oriNorm;
    double prev_error = relative_error+1;
    arma::mat update;
    int iterCount = 0, maxIter = 100;
    while(relative_error < prev_error && iterCount < maxIter){
        prev_error = relative_error;
        update=rInv*(-error);
        error= m_linkage.submat( start, start, endOfBlock, endOfBlock )*(heritResult+update) - fStat;
        relative_error = norm(error) / oriNorm;
        if(relative_error < 1e-300) relative_error = 0; // 1e-300 is more than enough...
        heritResult = heritResult+update;
        iterCount++; // This is to avoid infinite looping
    }
    arma::vec minusF = 1-fStat;
    for(size_t i = 0; i < fStat.n_elem; ++i){
        minusF(i) /= nSample(i)-2.0+zStat(i)*zStat(i);
    }
    // Now calculate the varaince, this is for the complicated variance
    varResult = rInv*arma::diagmat(minusF)*(4.0*m_linkageSqrt.submat(start,start,endOfBlock,endOfBlock)%(zStat*zStat.t()))*arma::diagmat(minusF)*rInv;
//    varResult = rInv*arma::diagmat(nSample)*(4.0*m_linkageSqrt.submat(start,start,endOfBlock,endOfBlock)%(zStat*zStat.t()))*arma::diagmat(nSample)*rInv;
    addVarResult = rInv*arma::diagmat(minusF)*(-2.0*m_linkage.submat(start,start,endOfBlock, endOfBlock))*arma::diagmat(minusF)*rInv;
//    addVarResult = rInv*arma::diagmat(nSample)*(-2.0*m_linkage.submat(start,start,endOfBlock, endOfBlock))*arma::diagmat(nSample)*rInv;
}


