#ifndef LINKAGE_H
#define LINKAGE_H


#include <algorithm>
#include <limits>
#include <math.h>


#include <boost/ptr_container/ptr_list.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <list>
#include <iterator>
#include <deque>

#include <stdio.h>
#include <iostream>
#include <assert.h>

#include <armadillo>

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
        void construct(boost::ptr_list<Genotype> &genotype, std::list<size_t> &snpLoc, std::deque<std::list<size_t>::iterator > &boundary, boost::ptr_vector<Snp> &snpList, const bool correction);
    protected:
    private:
        arma::mat m_linkage;
        arma::mat m_linkageSqrt;
        size_t m_thread=1;
        static std::mutex linkageMtx;
        // This will return the list of index that we would like to remove from the analysis
        void computeLd(const boost::ptr_list<Genotype> &genotype, const std::list<size_t> &snpLoc, size_t startIndex, size_t verEnd, size_t horistart,size_t horiEnd, boost::ptr_vector<Snp> &snpList, const bool &correction, std::vector<size_t> &perfectLd);
};

#endif // LINKAGE_H
