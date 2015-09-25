#include "decomposition.h"

//require major revision. There are BUGS here. So we need to make it robust to changes.

std::mutex Decomposition::mtx;
Decomposition::Decomposition(size_t thread):m_thread(thread){}

void Decomposition::solve(const std::vector<size_t> boundaries, const size_t index, const Linkage &linkage, const std::deque<size_t> snpLoc, boost::ptr_vector<Snp> &snpList, const bool &chromosomeStart, const Eigen::MatrixXd &betaEstimate){
    size_t sizeOfWindow = 0;
    size_t copyStart = boundaries[index+1]; //This is the index
    size_t copyEnd = copyStart; //This is the bound
    if(index+3==boundaries.size()){
        sizeOfWindow =linkage.size()-boundaries[index];
        copyEnd = linkage.size();
    }
    else{
        sizeOfWindow = boundaries[index+3]-boundaries[index]+1;
        copyEnd = boundaries[index+2];
    }
    if(chromosomeStart && index ==0){
        copyStart = boundaries[index];
    }
    Eigen::MatrixXd heritability = Eigen::MatrixXd::Zero(sizeOfWindow, 1);
    Eigen::MatrixXd effectiveNumber = Eigen::MatrixXd::Zero(sizeOfWindow, 1);
    Eigen::VectorXd ldScore = Eigen::VectorXd::Zero(sizeOfWindow);
    linkage.solve(boundaries[index], sizeOfWindow,betaEstimate, heritability,effectiveNumber,ldScore);
    //Now that we have got the results, we need to write them to the SNP
    for(size_t i = copyStart; i < copyEnd; ++i){
        snpList.at(snpLoc[i]).setHeritability(heritability(i-boundaries[index],0));
        snpList.at(snpLoc[i]).setEffectiveNumber(effectiveNumber(i-boundaries[index],0));
        snpList.at(snpLoc[i]).setLDScore(ldScore(i-copyStart));
    }
}

void Decomposition::run(const Linkage &linkage, const size_t &genotypeIndex, const size_t& remainedLD, const boost::ptr_vector<Interval> &blockInfo, std::deque<size_t> &ldLoc, const std::deque<size_t> &snpLoc, boost::ptr_vector<Snp> &snpList, const bool &chromosomeStart){
    size_t startRange =  genotypeIndex;
    size_t endRange=0;
    size_t i = genotypeIndex;
    size_t range = m_thread;
    if(remainedLD==0)range +=2;
    else{
        //Something was left behind, and they must be 2 blocks before
        startRange -= 2;//Two steps is only right if I have remove stuffs
    }
    std::string currentChr = blockInfo[startRange].getChr();
    for(; i< genotypeIndex+range && i < blockInfo.size(); ++i){
            if(blockInfo[i].getChr().compare(currentChr)!=0){
            i= i-1; //This mean we working on the last block of this chromosome
            break;
        }
    }
    endRange = i;
    std::vector<std::thread> threadStore;
    //Now every 3 blocks will be one window for decomposition
    //Launch a group of threads
    //Need start, end and mid bound
    //Now get all the boundaries
    std::vector<size_t> boundaries;
    size_t lastIndex =0;
    for(size_t i = startRange; i < endRange; ++i){
        for(; lastIndex < ldLoc.size(); ++lastIndex){
            if(ldLoc[lastIndex]==blockInfo[i].getStart()){
                boundaries.push_back(lastIndex);
                break;
            }
        }
    }
    Eigen::MatrixXd betaEstimate = Eigen::MatrixXd::Zero(ldLoc.size(),1);
    for(size_t i = 0; i < snpLoc.size(); ++i){
        betaEstimate(i)  = snpList.at(snpLoc[i]).getBeta();
    }
    //Now the boundaries should be representative of the location on the ld matrix
    size_t threadRunCounter =0;
    while(threadRunCounter < boundaries.size()-2){
        while(threadStore.size() < m_thread && threadRunCounter < boundaries.size()-2){ //On purposely leave 1 thread out for the main
            threadStore.push_back(std::thread(&Decomposition::solve, this, std::cref(boundaries), threadRunCounter, std::cref(linkage), std::cref(snpLoc), std::ref(snpList),chromosomeStart, std::cref(betaEstimate)));
            threadRunCounter++;
        }

        for (size_t j = 0; j < threadStore.size(); ++j) {
            threadStore[j].join();
        }
        threadStore.clear();
    }
}
