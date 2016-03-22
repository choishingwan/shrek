#include "decomposition.h"


Decomposition::Decomposition(size_t thread):m_thread(thread){}

Decomposition::~Decomposition()
{
    //dtor
}


void Decomposition::decompose(Linkage &linkage, std::deque<size_t> &snpLoc, size_t startDecompIter, size_t endDecompIter, size_t startVarIter, size_t endVarIter, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> &regionList, bool sign, bool start, bool ending){
    size_t sizeOfMatrix = endDecompIter-startDecompIter;
    arma::vec fStat = arma::vec(sizeOfMatrix, arma::fill::zeros);
    arma::vec heritResult = arma::vec(sizeOfMatrix, arma::fill::zeros);
    size_t regionSize = regionList.size();

    size_t startCopySnpLocIndex=(start)?startDecompIter : startVarIter;
    size_t startCopyVectorIndex = (start)? 0:startVarIter-startDecompIter;
    if(sign){
        arma::vec zStat = arma::vec(sizeOfMatrix, arma::fill::zeros);
        arma::vec nSample = arma::vec(sizeOfMatrix, arma::fill::zeros);
        arma::mat varResult = arma::mat(sizeOfMatrix,sizeOfMatrix,arma::fill::eye);
        size_t i = 0;
        for(size_t snpLocIndex = startDecompIter; snpLocIndex< endDecompIter; ++snpLocIndex){
            double stat = snpList[snpLoc[snpLocIndex]].getStat();
            int sampleSize = snpList[snpLoc[snpLocIndex]].getSampleSize();
            zStat(i) = stat;
            fStat(i) = snpList[snpLoc[snpLocIndex]].getFStat(); // This is the f-statistic
            nSample(i) = (double)1/(double)sampleSize;
//            nSample(i) = (double)sampleSize;
            ++i;
        }
        linkage.decompose(startDecompIter,zStat, fStat, nSample, heritResult, varResult);
        // Update the heritability (we separate it because it is easier)
        for(size_t j = startCopyVectorIndex; j < sizeOfMatrix; ++j){ // we don't care about the end in this case because they will be over wrote
            snpList[snpLoc[startCopySnpLocIndex]].setHeritability(heritResult(j));
            startCopySnpLocIndex++;
        }
        //The top block of the variance
        size_t horiCopyIndex = (start)? 0:endVarIter-startDecompIter;
        size_t startBound =startVarIter-startDecompIter;
        size_t endBound =endDecompIter-startDecompIter;
        for(size_t i =0; i < startBound; ++i){ // i is the index on the matrix
            for(size_t j = horiCopyIndex; j < endBound; ++j){ //  Again, j is the index on the matrix
                 for(size_t k = 0; k < regionSize; ++k){ // Iterate through the region list
                    if(snpList[snpLoc[startDecompIter+i]].flag(k) &&
                       snpList[snpLoc[startDecompIter+j]].flag(k))
                            regionList[k].addVariance(varResult(i,j));
                }
            }
        }
        // Now the mid block of the variance
        for(size_t i = startBound; i < endVarIter-startDecompIter; ++i){
            for(size_t j = 0; j < endBound; ++j){
                for(size_t k = 0; k < regionSize; ++k){ // Iterate through the region list
                    if(snpList[snpLoc[startDecompIter+i]].flag(k) &&
                       snpList[snpLoc[startDecompIter+j]].flag(k))
                            regionList[k].addVariance(varResult(i,j));
                }
            }
        }
        // Now the final block of the variance
        horiCopyIndex = (ending)? endBound: startBound;
        for(size_t i = endVarIter-startDecompIter; i < endBound; ++i){
            for(size_t j = 0; j < horiCopyIndex; ++j){
                for(size_t k = 0; k < regionSize; ++k){
                    if(snpList[snpLoc[startDecompIter+i]].flag(k) &&
                           snpList[snpLoc[startDecompIter+j]].flag(k))
                                regionList[k].addVariance(varResult(i,j));
                }
            }
        }
    }
    else{
        arma::vec varResult = arma::vec(sizeOfMatrix, arma::fill::zeros);
        size_t i = 0;
        for(size_t snpLocIndex = startDecompIter; snpLocIndex< endDecompIter; ++snpLocIndex){
            double stat = snpList[snpLoc[snpLocIndex]].getStat();
            fStat(i) = snpList[snpLoc[snpLocIndex]].getFStat();; // This is the f-statistic
            ++i;
        }
        linkage.decompose(startDecompIter,fStat,heritResult, varResult);
//        std::cerr << "Start copying: " << start << " " << startVarIter << " " << startDecompIter << std::endl;
        for(size_t j = startCopyVectorIndex; j < sizeOfMatrix; ++j){ // we don't care about the end in this case because they will be over wrote
            snpList[snpLoc[startCopySnpLocIndex]].setHeritability(heritResult(j));
            snpList[snpLoc[startCopySnpLocIndex]].setEffective(varResult(j));
            startCopySnpLocIndex++;
        }
//        std::cerr << "Finished copying" << std::endl;
    }
}


void Decomposition::run(Linkage &linkage, std::deque<size_t> &snpLoc, std::vector<size_t> &boundary, boost::ptr_vector<Snp> &snpList, bool windowEnd, bool decomposeAll, bool starting, boost::ptr_vector<Region> &regionList){
    // The proper decomposition process
    bool sign = !(snpList.front().getSign()==0); // Check whether if sign is given
    size_t boundSize = boundary.size();
    size_t snpLocSize = snpLoc.size();
    // Perform some basic checking
    if(!windowEnd) decomposeAll = false;
    assert(!(boundSize>4) && "Why? Not possible...");
    assert(!(starting && !windowEnd && boundSize != 4) && "Must have 4 boundaries if this is the start and not the end");
    assert(!(decomposeAll && !windowEnd) && "Decompose All only possible for the end!");
    assert(!(!starting && boundSize < 3+!windowEnd) &&"Insufficient block for decomposition");
    assert(!(!decomposeAll && windowEnd && boundSize==1)&& "Must decompose all if this is the end and only has one bound");

    if(boundSize ==1){
        // Check here
        assert((starting || windowEnd) && "How else is this possible?");
        decompose(linkage, snpLoc,
                  0, snpLocSize,
                  0, snpLocSize,
                  snpList, regionList, sign, true, true);
    }
    else if(boundSize ==2){
        assert((starting || windowEnd )&& "How else is this possible?");
        size_t midEnd = (decomposeAll)? snpLocSize:boundary.back();
        decompose(linkage, snpLoc,
                  0, midEnd,
                  0, midEnd,
                  snpList, regionList, sign, true, true);
    }
    else if(boundSize==3){
        assert((windowEnd) && "For this to happen, this must either be the end");

        if(!decomposeAll){
            assert(starting && "Must also be the start");
            decompose(linkage, snpLoc,
                  0, boundary[2],
                  boundary[1], boundary[2],
                  snpList, regionList, sign, starting, windowEnd);
        }
        else{
            decompose(linkage, snpLoc,
                  0, snpLocSize,
                  boundary[1], boundary[2],
                  snpList, regionList, sign, starting, windowEnd);
        }
    }
    else if(boundSize==4){
        decompose(linkage, snpLoc,
                  0, boundary.back(),
                  boundary[1], boundary[2],
                  snpList, regionList, sign, starting, windowEnd && !decomposeAll);
        if(decomposeAll && windowEnd){ //Only if this is windowEnd, otherwise, don't bother doing this
            decompose(linkage, snpLoc,
                      boundary[1], snpLocSize,
                      boundary[2], boundary.back(),
                      snpList, regionList, sign, starting, windowEnd && decomposeAll);
        }
    }
}
