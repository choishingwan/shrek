#include "decomposition.h"


Decomposition::Decomposition(size_t thread):m_thread(thread){}

Decomposition::~Decomposition()
{
    //dtor
}

void Decomposition::run(Linkage &linkage, std::list<size_t> &snpLoc, std::deque<std::list<size_t>::iterator > &boundary, boost::ptr_vector<Snp> &snpList, bool decomposeAll, size_t roundNumber, boost::ptr_vector<Region> &regionList){
    // Here is the meat of the decomposition preprocessing
    // The decomposition is actually done by the linkage class
    // First, get the vectors ready
    arma::vec fStat = arma::vec(snpLoc.size(), arma::fill::zeros);
    arma::vec zStat = arma::vec(snpLoc.size(), arma::fill::zeros);
    arma::vec nSample = arma::vec(snpLoc.size(), arma::fill::zeros);
    std::list<size_t>::iterator snpLocIter = snpLoc.begin();
    size_t i = 0;
    bool checked = false, sign=false;
    for(; i < snpLoc.size(); ++i){
        double stat = snpList.at(*(snpLocIter)).getStat();
        int sampleSize = snpList.at(*(snpLocIter)).getSampleSize();
        if(!checked){
            sign = (snpList.at(*(snpLocIter)).getSign()==0)? false: true;
        }
        zStat(i) = stat;
        fStat(i) = (stat*stat-1.0)/((double)sampleSize-2.0+stat*stat); // This is the f-statistic
        nSample(i) = sampleSize;
        std::advance(snpLocIter, 1);
        if(!decomposeAll && snpLocIter == boundary.back()) break;
    }
    if(!decomposeAll){
        zStat.resize(i);
        fStat.resize(i);
        nSample.resize(i);
    }
    /** Now we have all the required factors **/
    arma::vec heritResult = arma::vec(i, arma::fill::zeros);

    // These are for updating stuffs
    std::list<size_t>::iterator iter = snpLoc.begin();
    size_t numItem = i; //For now, i still represent the size of the result vector
    if(roundNumber == 0) i = 0; // this is the start of the region, we will want to get all the results
    else{
        i = std::distance(snpLoc.begin(), boundary[1]); // It means we will start at the middle
        iter = boundary[1];
    }

    if(sign){
        // Because we have the sign, we will perform the more complicated variance calculation
        arma::mat varResult = arma::mat(i,i,arma::fill::eye);
        linkage.decompose(zStat, fStat, nSample, heritResult, varResult);
        // Special variance handling here
    }
    else{
        arma::vec varResult = arma::vec(i, arma::fill::zeros);
        linkage.decompose(fStat,heritResult, varResult);
        // In this case, we can just update the effective number
        for(size_t j = i; j < numItem; ++j){
            snpList.at(*(iter)).setEffective(varResult(j));
            std::advance(iter, 1);
        }
    }
    // The update of heritability will always be the same for both


    // As now we don't perform multi-threading here, it is safe to just update all the ending stuff without worrying
    // about out of bound stuff
    if(roundNumber != 0) iter = boundary[1]; // reset the iterator
    for(size_t j = i; j < numItem; ++j){
        snpList.at(*(iter)).setHeritability(heritResult(j));
        std::advance(iter, 1);
    }

}
