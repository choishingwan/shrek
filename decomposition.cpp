#include "decomposition.h"


Decomposition::Decomposition(size_t thread):m_thread(thread){}

Decomposition::~Decomposition()
{
    //dtor
}


void Decomposition::decompose(Linkage &linkage, std::deque<size_t> &snpLoc, size_t startDecompIter, size_t endDecompIter, size_t startVarIter, size_t endVarIter, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> &regionList, bool sign, bool start, bool ending){
    size_t sizeOfMatrix = endDecompIter-startDecompIter;
    size_t startCoordinate = startDecompIter;
    arma::vec fStat = arma::vec(sizeOfMatrix, arma::fill::zeros);
    arma::vec heritResult = arma::vec(sizeOfMatrix, arma::fill::zeros);
    size_t regionSize = regionList.size();
    if(sign){
        arma::vec zStat = arma::vec(sizeOfMatrix, arma::fill::zeros);
        arma::vec nSample = arma::vec(sizeOfMatrix, arma::fill::zeros);
        arma::mat varResult = arma::mat(sizeOfMatrix,sizeOfMatrix,arma::fill::eye);
        arma::mat addVarResult = arma::mat(sizeOfMatrix,sizeOfMatrix,arma::fill::eye);
        size_t i = 0;
        for(size_t snpLocIndex = startDecompIter; snpLocIndex< endDecompIter; ++snpLocIndex){
            double stat = snpList.at(snpLoc[snpLocIndex]).getStat();
            int sampleSize = snpList.at(snpLoc[snpLocIndex]).getSampleSize();
            zStat(i) = stat;
            fStat(i) = (stat*stat-1.0)/((double)sampleSize-2.0+stat*stat); // This is the f-statistic
            nSample(i) = (double)1/(double)sampleSize;
//            nSample(i) = (double)sampleSize;
            ++i;
        }
        linkage.decompose(startDecompIter,zStat, fStat, nSample, heritResult, varResult);
//        linkage.decompose(startCoordinate,zStat, fStat, nSample, heritResult, varResult, addVarResult);
        // Update the heritability (we separate it because it is easier)
        size_t startCopySnpLocIndex=(start)?startDecompIter : startVarIter;
        size_t startCopyVectorIndex = (start)? 0:startVarIter-startDecompIter;
        for(size_t j = startCopyVectorIndex; j < sizeOfMatrix; ++j){ // we don't care about the end in this case because they will be over wrote
            snpList.at(snpLoc[startCopySnpLocIndex]).setHeritability(heritResult(j));
            startCopySnpLocIndex++;
        }
        //The top block of the variance
        size_t horiCopyIndex = (start)? 0:endVarIter-startDecompIter;
        size_t startBound =startVarIter-startDecompIter;
        size_t endBound =endDecompIter-startDecompIter;
        for(size_t i =0; i < startBound; ++i){ // i is the index on the matrix
            for(size_t j = horiCopyIndex; j < endBound; ++j){ //  Again, j is the index on the matrix
                 for(size_t k = 0; k < regionSize; ++k){ // Iterate through the region list
                    if(snpList.at(snpLoc[startDecompIter+i]).flag(k) &&
                       snpList.at(snpLoc[startDecompIter+j]).flag(k))
                            regionList[k].addVariance(varResult(i,j));
//                            regionList[k].addVariance(varResult(i,j), addVarResult(i,j));
                }
            }
        }
        // Now the mid block of the variance
        for(size_t i = startBound; i < endVarIter-startDecompIter; ++i){
            for(size_t j = 0; j < endBound; ++j){
                for(size_t k = 0; k < regionSize; ++k){ // Iterate through the region list
                    if(snpList.at(snpLoc[startDecompIter+i]).flag(k) &&
                       snpList.at(snpLoc[startDecompIter+j]).flag(k))
                            regionList[k].addVariance(varResult(i,j));
//                            regionList[k].addVariance(varResult(i,j), addVarResult(i,j));
                }
            }
        }
        // Now the final block of the variance
        horiCopyIndex = (ending)? endBound: startBound;
        for(size_t i = endVarIter-startDecompIter; i < endBound; ++i){
            for(size_t j = 0; j < horiCopyIndex; ++j){
                for(size_t k = 0; k < regionSize; ++k){
                    if(snpList.at(snpLoc[startDecompIter+i]).flag(k) &&
                           snpList.at(snpLoc[startDecompIter+j]).flag(k))
                                regionList[k].addVariance(varResult(i,j));
//                                regionList[k].addVariance(varResult(i,j), addVarResult(i,j));
                }
            }
        }
    }
    else{
        arma::vec varResult = arma::vec(sizeOfMatrix, arma::fill::zeros);
        size_t i = 0;
        for(size_t snpLocIndex = startDecompIter; snpLocIndex< endDecompIter; ++snpLocIndex){
            double stat = snpList.at(snpLoc[snpLocIndex]).getStat();
            int sampleSize = snpList.at(snpLoc[snpLocIndex]).getSampleSize();
            fStat(i) = (stat*stat-1.0)/((double)sampleSize-2.0+stat*stat); // This is the f-statistic
            ++i;
        }
        linkage.decompose(startDecompIter,fStat,heritResult, varResult);
        size_t startCopySnpLocIndex=(start)?startDecompIter : startVarIter;
        size_t startCopyVectorIndex = (start)? 0:startVarIter-startDecompIter;
        for(size_t j = startCopyVectorIndex; j < sizeOfMatrix; ++j){ // we don't care about the end in this case because they will be over wrote
            snpList.at(snpLoc[startCopySnpLocIndex]).setHeritability(heritResult(j));
            snpList.at(snpLoc[startCopySnpLocIndex]).setEffective(varResult(j));
            startCopySnpLocIndex++;
        }
    }
}


void Decomposition::run(Linkage &linkage, std::deque<size_t> &snpLoc, std::deque<size_t> &boundary, boost::ptr_vector<Snp> &snpList, bool finalizeBuff, bool decomposeAll, bool starting, boost::ptr_vector<Region> &regionList){
    // The proper decomposition process
    bool sign = !(snpList.front().getSign()==0); // Check whether if sign is given
    size_t boundSize = boundary.size();
    size_t snpLocSize = snpLoc.size();
    // Perform some basic checking

    assert(!(boundSize>4) && "Why? Not possible...");
    assert(!(starting && !finalizeBuff && boundSize != 4) && "Must have 4 boundaries if this is the start and not the end");
    assert(!(decomposeAll && !finalizeBuff) && "DecomposeAll only possible for the end!");
    assert(!(!starting && boundSize < 3+!finalizeBuff) &&"Insufficient block for decomposition");
    assert(!(!decomposeAll && finalizeBuff && boundSize==1)&& "Must decompose all if this is the end and only has one bound");

    if(boundSize ==1){
        // Check here
        assert((starting || finalizeBuff) && "How else is this possible?");
        decompose(linkage, snpLoc,
                  0, snpLocSize,
                  0, snpLocSize,
                  snpList, regionList, sign, true, true);
    }
    else if(boundSize ==2){
        assert((starting || finalizeBuff )&& "How else is this possible?");
        size_t midEnd = (decomposeAll)? snpLocSize:boundary.back();
        decompose(linkage, snpLoc,
                  0, midEnd,
                  0, midEnd,
                  snpList, regionList, sign, true, true);
    }
    else if(boundSize==3){
        assert((starting||finalizeBuff) && "For this to happen, this must either be the starting or the end");
        if(!decomposeAll){
            assert(starting && finalizeBuff && "Must be the start and the end!");
            decompose(linkage, snpLoc,
                  0, boundary[2],
                  0, boundary[2],
                  snpList, regionList, sign, starting, finalizeBuff);
        }
        else{
            assert(finalizeBuff && "Must be the end!");
            decompose(linkage, snpLoc,
                  0, snpLocSize,
                  0, boundary[2],
                  snpList, regionList, sign, starting, finalizeBuff);
        }
    }
    else if(boundSize==4){
        decompose(linkage, snpLoc,
                  0, boundary.back(),
                  boundary[1], boundary[2],
                  snpList, regionList, sign, starting, finalizeBuff && !decomposeAll);
        if(decomposeAll && finalizeBuff){ //Only if this is finalizeBuff, otherwise, don't bother doing this
            decompose(linkage, snpLoc,
                      boundary[1], snpLocSize,
                      boundary[2], boundary.back(),
                      snpList, regionList, sign, starting, finalizeBuff && decomposeAll);
        }
    }
}
