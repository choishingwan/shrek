#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include <iterator>
#include <deque>
#include <list>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <armadillo>

#include "linkage.h"
#include "snp.h"

class Decomposition
{
    public:
        /** Default constructor */
        Decomposition(size_t thread);
        /** Default destructor */
        virtual ~Decomposition();
        void run(Linkage &linkage, std::list<size_t> &snpLoc, std::deque<std::list<size_t>::iterator > &boundary, boost::ptr_vector<Snp> &snpList, bool decomposeAll, size_t roundNumber, boost::ptr_vector<Region> &regionList);
    protected:
    private:
        size_t m_thread = 0;

};

#endif // DECOMPOSITION_H
