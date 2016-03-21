#ifndef LINKAGE_H
#define LINKAGE_H


#include <algorithm>
#include <limits>
#include <math.h>

#include <stdexcept>

//#include <boost/ptr_container/ptr_list.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <list>
#include <iterator>
#include <deque>
#include <algorithm>

#include <stdio.h>
#include <iostream>
#include <assert.h>

#include <armadillo>
#include <Eigen/Dense>
#include <limits>
#include <mutex>
#include <thread>

#include "genotype.h"
#include "snp.h"

class Linkage
{
    public:
        /** Default constructor */
		Linkage(size_t thread);
        /** Default destructor */
        virtual ~Linkage();
        // The snpList is required for the perfectLD stuff
        void construct(boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::vector<size_t> &boundary, boost::ptr_vector<Snp> &snpList, const bool correction,bool &boundCheck);
        void print();
        void decompose(size_t start, const arma::vec &fStat, arma::vec &heritResult, arma::vec &varResult);
        void decompose(size_t start, const arma::vec &zStat, const arma::vec &fStat, const arma::vec &nSample, arma::vec &heritResult, arma::mat &varResult);
        void computeHerit(const arma::vec &fStat, arma::vec &heritResult);
        void clear();
        void clear(size_t nRemoveElements);
    protected:
    private:
        arma::mat m_linkage;
        arma::mat m_linkageSqrt;
        size_t m_thread=1;
        size_t m_blockSize=0;
        static std::mutex linkageMtx;
        // This will return the list of index that we would like to remove from the analysis
        void computeLd(const boost::ptr_deque<Genotype> &genotype, const std::deque<size_t> &snpLoc, size_t startIndex, size_t verEnd, size_t horistart, boost::ptr_vector<Snp> &snpList, const bool &correction, std::vector<size_t> &perfectLd);
        void perfectRemove(std::vector<size_t> &perfectLd, boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::vector<size_t > &boundary, boost::ptr_vector<Snp> &snpList, bool &boundCheck);
};

#endif // LINKAGE_H
