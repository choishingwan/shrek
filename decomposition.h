#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include <vector>
#include "snpindex.h"
#include "snp.h"
#include "linkage.h"

class Decomposition
{
	public:
		Decomposition(SnpIndex *snpIndex, std::vector<Snp*> *snpList, Linkage *linkageMatrix, size_t thread);
		virtual ~Decomposition();
	protected:
	private:
		SnpIndex *m_snpIndex;
        std::vector<Snp*> *m_snpList;
        Linkage *m_linkage;
		size_t m_thread;

};

#endif // DECOMPOSITION_H
