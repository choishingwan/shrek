#include "linkage.h"

std::mutex Linkage::linkageMtx;

Linkage::Linkage(size_t thread, size_t blockSize):m_thread(thread), m_blockSize(blockSize){}
Linkage::~Linkage(){}

void Linkage::print(){ std::cout << m_linkage << std::endl; }

void Linkage::computeLd(const boost::ptr_list<Genotype> &genotype, const std::list<size_t> &snpLoc, size_t startIndex, size_t verEnd, size_t horistart, size_t horiEnd, boost::ptr_vector<Snp> &snpList, const bool &correction, std::vector<size_t> &perfectLd){
    std::vector<size_t> perfectBuff;
    boost::ptr_list<Genotype>::const_iterator  iter = genotype.begin(), endIter = genotype.begin(), horiEndIter=genotype.begin();
    std::advance(iter, startIndex);
    std::advance(endIter, verEnd);
    std::advance(horiEndIter, horiEnd);

    size_t i = startIndex;
    for(; iter != endIter; ++iter){
        size_t nextStart = (i>horistart)? i:horistart;
        size_t j = nextStart;
        boost::ptr_list<Genotype>::const_iterator  horiStartIter=genotype.begin();
        std::advance(horiStartIter, nextStart);
        for(; horiStartIter!=horiEndIter;  ++horiStartIter){
                if(m_linkage(i,j)==0.0){  // only calculate the R2 when this is a new spot
                    if(i==j){
                        m_linkage(i,j) = 1.0;
                        m_linkageSqrt(i,j) = 1.0;
                    }
                    else{
                        double r = 0.0;
                        double r2 = 0.0;
                        (*iter).GetbothR((*horiStartIter), correction, r, r2);
                        m_linkage(i,j) = r2;
                        m_linkageSqrt(i,j) = r;
                        // Now check if it is perfect LD
                        // Store the index j into perfectBuff
                        if(std::fabs(r-1.0) < 1e-6 || r > 1.0 || r < -1.0){
                            perfectBuff.push_back(j);
                            //This can be slow?? But then, this should also be rare
                            std::list<size_t>::const_iterator currentSnp = snpLoc.begin();
                            std::advance(currentSnp, std::distance(genotype.begin(), iter));
                            std::list<size_t>::const_iterator targetSnp = snpLoc.begin();
                            std::advance(targetSnp, std::distance(genotype.begin(), horiStartIter));
                            snpList.at((*currentSnp)).shareHeritability(snpList.at((*targetSnp)));
                        }
                    }
                }
            //This is for the update of the matrix indexing;
            ++j;
        }
        ++i; //Again, this is for the update of matrix index
    }
    // Now we can insert the buff back to the result
    linkageMtx.lock();
        perfectLd.insert(perfectLd.end(), perfectBuff.begin(), perfectBuff.end());
    linkageMtx.unlock();

}


void Linkage::construct(boost::ptr_list<Genotype> &genotype, std::list<size_t> &snpLoc, std::deque<std::list<size_t>::iterator > &boundary, boost::ptr_vector<Snp> &snpList, const bool correction,bool &boundCheck){
    // This is yet another complicated algorithm
    // The aim of this function is to construct the required LD matrix
    // (Think about it, maybe I will shift the perfect LD removal to the other function)
    if(genotype.empty()){
        throw std::runtime_error("Cannot build LD without genotypes");
	}
    // If this is the first of its kind, then we need to initialize the matrix
//    std::cerr << "Checking: " << m_linkage.n_cols << "\t" << genotype.size() << std::endl;
	if(m_linkage.n_cols==0){
        m_linkage = arma::mat(genotype.size(), genotype.size(),arma::fill::eye);
        m_linkageSqrt = arma::mat(genotype.size(), genotype.size(), arma::fill::eye);
    }
    else{ // Otherwise, just make sure the matrix is of the right size
        m_linkage.resize(genotype.size(), genotype.size());
        m_linkageSqrt.resize(genotype.size(), genotype.size());
    }
//    std::cerr << "Size expanded " << boundary.size() << std::endl;
    // We don't need the buff information
    // because for this function, we should build all the matrix anyway
    // Now if we are not working on the new stuff, we assume the LD matrix has been clean
    // of the unwanted materials
    // So we only need to construct the matrix here
    // When we construct the LD matrix, we do it per SNP level but still use multi-threading
    // then we skip any SNP that are found to be in perfect LD with others
    // This will allow us to avoid the need of removing rols (only column)
    // which should be fast considering armadillo is column major
    std::vector<size_t> perfectLd; // This is the vector use to store the perfect LD information
    std::vector<std::thread> threadStore;
    if(boundary.size() < 4){
        //Very small input, so we will do the whole part together
        size_t startIndex = 0, endBound= genotype.size();
        if(genotype.size() < m_thread){
            // Use as much thread as we can
            for(size_t i = startIndex; i < endBound; ++i){
                threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), i, i+1, i,genotype.size(), std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
            }
            for (size_t i = 0; i < threadStore.size(); ++i) {
                threadStore[i].join();
            }
            threadStore.clear();
        }
        else{
            // distribute the work
            int jobBlock = endBound/m_thread;
            int remainder = endBound % m_thread;
            int currentBlock = startIndex;
            for(size_t i = 0; i < m_thread; ++i){
                // For each thread
                if(remainder > 0){
                    threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), currentBlock, currentBlock+jobBlock+1, currentBlock, genotype.size(), std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
                    --remainder;
                    currentBlock+=jobBlock+1;
                }
                else{
                    threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), currentBlock, currentBlock+jobBlock, currentBlock, genotype.size(), std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
                    currentBlock+=jobBlock;
                }

            }
            for (size_t i = 0; i < threadStore.size(); ++i) {
                threadStore[i].join();
            }
            threadStore.clear();
        }
    }
    else{ // This is big enough, more importantly, not we have the extra bit for us to work with
          // With the new method, it will never be start
          // Therlinkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);efore, the first block are always completed
//        std::cerr << "Trying to build the last block" << std::endl;



        // Important thing, startIndex, verEnd, horiStart and horiEnd are all index of snpList and are not corresponding to the LD matrix location
        // rather, their distance from the beginning of snpLoc is


        size_t startIndex = std::distance(snpLoc.begin(), boundary[1]); // Start from the second block
        size_t endBound = genotype.size();
        // We do have enough stuff to work on
        int numSnps = endBound-startIndex;
        assert(numSnps>0 && "There must be some SNPs for the construct to work with!");

        size_t horiStart = std::distance(snpLoc.begin(), boundary.back());
        if(numSnps< (signed) m_thread){
            //very unlikely, but we will consider it anyway
            for(size_t i = startIndex; i < endBound; ++i){
//                std::cerr << i << " " << i+1 << " " << horiStart << " " << endBound << std::endl;
                threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), i, i+1, horiStart, endBound, std::ref(snpList), std::cref(correction), std::ref(perfectLd)));

            }
            for (size_t i = 0; i < threadStore.size(); ++i) {
                threadStore[i].join();
            }
            threadStore.clear();
        }
        else{
            // Divide the jobs
            int jobBlock = numSnps/m_thread;
            int remainder =numSnps%m_thread;
            int currentBlock = startIndex;
            for(size_t i = 0; i < m_thread; ++i){
                if(remainder>0){
                    threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), currentBlock, currentBlock+jobBlock+1, horiStart, endBound, std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
//                    std::cerr << currentBlock << " " << currentBlock+jobBlock+1 << " " << horiStart << " " << endBound << std::endl;
                    currentBlock+= jobBlock+1;
                    remainder--;
                }
                else{
                    threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), currentBlock, currentBlock+jobBlock, horiStart, endBound, std::ref(snpList), std::cref(correction), std::ref(perfectLd)));
//                    std::cerr << currentBlock << " " << currentBlock+jobBlock << " " << horiStart << " " << endBound << std::endl;
                    currentBlock+= jobBlock;
                }
            }
            for (size_t i = 0; i < threadStore.size(); ++i) {
                threadStore[i].join();
            }
            threadStore.clear();
        }
        // Join all the results

    }

    // Now we have the information of where all the perfect LD locates
    std::sort(perfectLd.begin(), perfectLd.end());
    // Now we want to remove duplicated index
    perfectLd.erase( std::unique( perfectLd.begin(), perfectLd.end() ), perfectLd.end() );
    // Here goes the crazy part
//     Need to remove the perfect LD
    // Now we need to handle the perfect LD stuff
    if(perfectLd.size()==0){
        // there is nothing to work on
        m_linkage = arma::symmatu(m_linkage);
        m_linkageSqrt = arma::symmatu(m_linkageSqrt);
    }
    else perfectRemove(perfectLd, genotype, snpLoc, boundary, snpList, boundCheck);
}


void Linkage::perfectRemove(std::vector<size_t> &perfectLd, boost::ptr_list<Genotype> &genotype, std::list<size_t> &snpLoc, std::deque<std::list<size_t>::iterator > &boundary, boost::ptr_vector<Snp> &snpList,bool &boundCheck){
// Challenge here: If the boundary is the perfect LD
// then we need to reconstruct stuff, also, if the next
// boundary is too far, we will have to cut the block...
    // No matter what, we need to update the SNP test statistics
    std::cerr << "m_linkage: " << std::endl;
    for(std::list<size_t>::iterator iter = snpLoc.begin(); iter != snpLoc.end(); ++iter) std::cerr << *(iter) << " ";
    std::cerr << std::endl;

    std::cerr << m_linkageSqrt << std::endl;
    for(size_t i = 0; i < boundary.size(); ++i) std::cerr << "Bound check: "<< *(boundary[i]) << std::endl;
    if(boundary.size() !=1){
        // we have to check the boundary, but only with the last one
        // perfectLD is sorted in ascending order
        for(size_t i = 0; i < perfectLd.size(); ++i){
            std::list<size_t>::iterator currentSnp = snpLoc.begin();
            std::advance(currentSnp, perfectLd[i]);
            std::cerr <<"Checking: " << perfectLd[i] << std::endl;
            // We need to check, if the current bound is the last bound, then we
            // should just delete it?
            if(currentSnp == boundary.back()){
                // The boundary SNP is also in perfect LD with previous SNPs
                // So we will advance the pointer.
                std::advance(boundary.back(),1); // we will move the SNP to the next one that is not perfect LD
                // However, if the boundary is now at the end of the snpLoc, it means there is no SNP in this block that
                // is not in perfect LD
                // in this case, we will just remove the boundary?
                boundCheck = true;
                if(boundary.back()==snpLoc.end()){
                    break;
                }
            }
            // We will loop through the whole thing, such that we will make sure the boundary.back is pointing to
            // somewhere that is not in perfect LD with any existing SNPs
        }
        // We don't care whether if the blocks are alright, we just need to remove all the perfect LD here
    }
    // When we reach here, the perfect Ld stuff is solved and we should be able to just update the whole matrix
    // now update the matrix


    // Need to update a few things,
    // 1. Genotype
    // 2. snpLoc
    // 3. m_linkage + m_linageSqrt
    // Do it slowly here first, optimize later when I am certain this works
    size_t cI = 0;
    size_t perfectIndexI = 0;
    for(size_t i = 0; i < genotype.size(); ++i){
        if(i == perfectLd[perfectIndexI]) perfectIndexI++;// This is perfect LD, we will skip it
        else{
            // the current I is required, so we will go over the j

            size_t cJ = cI, perfectIndexJ = perfectIndexI;
            for(size_t j = i; j < genotype.size(); ++j){
                if(j == perfectLd[perfectIndexJ]) perfectIndexJ++;
                else{
                    m_linkage(cJ,cI) = m_linkage(i,j);
                    m_linkageSqrt(cJ,cI) = m_linkageSqrt(i,j);
                    cJ++;
                }
            }
            cI++;
        }
    }
    // now update the genotype and snpLoc
    std::vector<std::list<size_t>::iterator> snpLocIterRemove;
    std::vector<boost::ptr_list<Genotype>::iterator> genoIterRemove;
    for(size_t i = 0; i < perfectLd.size(); ++i){
        std::list<size_t>::iterator currentLoc = snpLoc.begin();
        boost::ptr_list<Genotype>::iterator currentGeno = genotype.begin();
        std::advance(currentLoc, perfectLd[i]);
        std::advance(currentGeno, perfectLd[i]);
        snpLocIterRemove.push_back(currentLoc);
        genoIterRemove.push_back(currentGeno);
    }
    for(size_t i = 0; i < perfectLd.size(); ++i){
        snpLoc.erase(snpLocIterRemove[i]);
        genotype.erase(genoIterRemove[i]);
    }
    //The vectors will go out of scope and they are not pointers, so don't bother cleaning them up
    m_linkage.resize(genotype.size(), genotype.size());
    m_linkageSqrt.resize(genotype.size(), genotype.size());
    m_linkage = arma::symmatl(m_linkage);
    m_linkageSqrt = arma::symmatl(m_linkageSqrt);

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
    size_t endOfBlock = fStat.n_elem-1;
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
}



