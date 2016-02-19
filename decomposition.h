#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include <iterator>
#include <deque>
#include <list>
#include <stdexcept>
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
//        void run(Linkage &linkage, std::list<size_t> &snpLoc, std::deque<std::list<size_t>::iterator > &boundary, boost::ptr_vector<Snp> &snpList, bool finalizeBuff, bool decomposeAll, bool starting, boost::ptr_vector<Region> &regionList);
        void run(Linkage &linkage, std::deque<size_t> &snpLoc, std::deque<size_t> &boundary, boost::ptr_vector<Snp> &snpList, bool finalizeBuff, bool decomposeAll, bool starting, boost::ptr_vector<Region> &regionList);
    protected:
    private:
        size_t m_thread = 0;
//        void decompose(Linkage &linkage, std::list<size_t> &snpLoc, std::list<size_t>::iterator startDecompIter, std::list<size_t>::iterator endDecompIter, std::list<size_t>::iterator startCopyIter, std::list<size_t>::iterator endCopyIter, std::list<size_t>::iterator startVarIter, std::list<size_t>::iterator endVarIter, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> &regionList, bool sign, bool start, bool ending);
        void decompose(Linkage &linkage, std::deque<size_t> &snpLoc, size_t startDecompIter, size_t endDecompIter, size_t startCopyIter, size_t endCopyIter, size_t startVarIter, size_t endVarIter, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> &regionList, bool sign, bool start, bool ending);

};

#endif // DECOMPOSITION_H
