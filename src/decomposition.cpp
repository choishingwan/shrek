#include "decomposition.h"


Decomposition::Decomposition(size_t thread):m_thread(thread){}

Decomposition::~Decomposition()
{
    //dtor
}


void Decomposition::complexSEUpdate(size_t decompStartIndex, size_t decompEndIndex, size_t midStartIndex, size_t midEndIndex, std::deque<size_t> &snpLoc,boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> &regionList, bool starting, bool ending, const arma::mat &varResult){
    size_t horiCopyIndex = (starting)? 0:midEndIndex-decompStartIndex;
    size_t startBound =midStartIndex-decompStartIndex;
    size_t endBound =decompEndIndex-decompStartIndex;
    size_t regionSize = regionList.size();
    for(size_t i =0; i < startBound; ++i){ // i is the index on the matrix
        for(size_t j = horiCopyIndex; j < endBound; ++j){ //  Again, j is the index on the matrix
            for(size_t k = 0; k < regionSize; ++k){ // Iterate through the region list
                if(snpList[snpLoc[decompStartIndex+i]].flag(k) &&
                        snpList[snpLoc[decompStartIndex+j]].flag(k))
                    regionList[k].addVariance(varResult(i,j));
            }
        }
    }
        // Now the mid block of the variance
    for(size_t i = startBound; i < midEndIndex-decompStartIndex; ++i){
        for(size_t j = 0; j < endBound; ++j){
            for(size_t k = 0; k < regionSize; ++k){ // Iterate through the region list
                if(snpList[snpLoc[decompStartIndex+i]].flag(k) &&
                    snpList[snpLoc[decompStartIndex+j]].flag(k))
                        regionList[k].addVariance(varResult(i,j));
            }
        }
    }
        // Now the final block of the variance
    horiCopyIndex = (ending)? endBound: startBound;
    for(size_t i= midEndIndex-decompStartIndex; i < endBound; ++i){
        for(size_t j = 0; j < horiCopyIndex; ++j){
            for(size_t k = 0; k < regionSize; ++k){
                if(snpList[snpLoc[decompStartIndex+i]].flag(k) &&
                        snpList[snpLoc[decompStartIndex+j]].flag(k))
                            regionList[k].addVariance(varResult(i,j));
            }
        }
    }
}

size_t Decomposition::check = 0;
void Decomposition::decompose(Linkage &linkage, std::deque<size_t> &snpLoc, const size_t decompStartIndex, const size_t decompEndIndex, const size_t midStartIndex, const size_t midEndIndex, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> &regionList, bool sign, bool starting, bool ending){
    size_t sizeOfMatrix = decompEndIndex-decompStartIndex;
    arma::vec fStat = arma::vec(sizeOfMatrix, arma::fill::zeros);
    arma::vec heritResult = arma::vec(sizeOfMatrix, arma::fill::zeros);
    arma::vec tStat, nSample;
    if(sign){
        tStat = arma::vec(sizeOfMatrix, arma::fill::zeros);
        nSample = arma::vec(sizeOfMatrix, arma::fill::zeros);
    }
    size_t regionSize = regionList.size();
    size_t startCopySnpLocIndex=(starting)? decompStartIndex : midStartIndex;
    size_t startCopyVectorIndex = (starting)? 0:midStartIndex-decompStartIndex;
    for(size_t snpLocIndex = decompStartIndex; snpLocIndex< decompEndIndex; ++snpLocIndex){
        fStat(snpLocIndex-decompStartIndex) = snpList[snpLoc[snpLocIndex]].getFStat(); // This is the f-statistic
        if(sign){
            tStat(snpLocIndex-decompStartIndex) = snpList[snpLoc[snpLocIndex]].getTStat();
            nSample(snpLocIndex-decompStartIndex) = (double)1/(double)snpList[snpLoc[snpLocIndex]].getSampleSize();
        }
    }

    linkage.decompose(decompStartIndex, decompEndIndex, fStat, heritResult);
    for(size_t i = (starting)? 0:midStartIndex-decompStartIndex; i < sizeOfMatrix; ++i){
        snpList[snpLoc[i+decompStartIndex]].setHeritability(heritResult(i));
    }
    if(sign){
        arma::mat varResult = arma::mat(sizeOfMatrix,sizeOfMatrix,arma::fill::eye);
        linkage.complexSE(decompStartIndex, decompEndIndex, nSample, tStat, varResult);
        complexSEUpdate(decompStartIndex, decompEndIndex, midStartIndex, midEndIndex, snpLoc, snpList, regionList, starting, ending, varResult);
    }
    else{
        arma::vec varResult = arma::vec(sizeOfMatrix, arma::fill::zeros);
        linkage.effectiveSE(decompStartIndex, decompEndIndex, varResult);
        for(size_t i = (starting)? 0:midStartIndex-decompStartIndex; i < sizeOfMatrix; ++i){
            snpList[snpLoc[i+decompStartIndex]].setEffective(heritResult(i));
        }
    }

}


void Decomposition::run(Linkage &linkage, std::deque<size_t> &snpLoc, std::vector<size_t> &boundary, boost::ptr_vector<Snp> &snpList, bool windowEnd, bool decomposeAll, bool starting, boost::ptr_vector<Region> &regionList){
    // The proper decomposition process
    bool sign = !(snpList.front().getSign()==0); // Check whether if sign is given
    size_t boundSize = boundary.size();
    size_t snpLocSize = snpLoc.size();
    // Perform some basic checking
    if(!windowEnd) decomposeAll = false;
    assert(!(boundSize>5) && "Why? Not possible...");
    assert(!(starting && !windowEnd && boundSize != 5) && "Must have 5 boundaries if this is the start and not the end");


    size_t decomposeStart = 0;
    size_t decomposeEnd = (decomposeAll)? ((boundSize>3)?boundary[3]:snpLocSize) : ((boundSize>3)? boundary[3]:boundary.back());
    assert(decomposeEnd!=decomposeStart && "Not possible to have 1 block and still want to retain stuff");
    size_t copyStart = (boundSize>3)? boundary[1]: boundary.front();
    size_t copyEnd = (boundSize>3)? boundary[2]: decomposeEnd;

    decompose(linkage, snpLoc,
            decomposeStart, decomposeEnd,
            copyStart, copyEnd,
            snpList, regionList, sign, starting,
            (boundSize>3&&(!decomposeAll&&windowEnd ))||(!boundSize>3 && windowEnd) );
    //Now handle the remaining
    if(windowEnd && boundSize>3){
        //Maximum two more situation
        for(size_t i = 0; i < boundSize-3-!decomposeAll;++i){
            decomposeStart=boundary[i+1];
            decomposeEnd = (i+4>=boundSize)? snpLocSize:boundary[i+4];
            copyStart = boundary[i+2];
            copyEnd = boundary[i+3];

            decompose(linkage, snpLoc,
                    decomposeStart, decomposeEnd,
                    copyStart, copyEnd,
                    snpList, regionList, sign, false,
                    ((i+1)==(boundSize-3-!decomposeAll)));
        }

    }
}
