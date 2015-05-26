#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include <vector>
#include <deque>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "snpindex.h"
#include "snp.h"
#include "region.h"
#include "linkage.h"
#include "genotype.h"
#include "processcode.h"
#include "decompositionthread.h"

class Decomposition
{
	public:
		Decomposition(SnpIndex *snpIndex, std::vector<Snp*> *snpList, Linkage *linkageMatrix, size_t thread);
		virtual ~Decomposition();
		ProcessCode Decompose(const size_t &blockSize, std::deque<size_t> &snpLoc, std::deque<Genotype*> &genotype, bool chromosomeStart, bool chromosomeEnd);
	protected:
	private:
		SnpIndex *m_snpIndex;
        std::vector<Snp*> *m_snpList;
        Linkage *m_linkage;
		size_t m_thread;

};

#endif // DECOMPOSITION_H
