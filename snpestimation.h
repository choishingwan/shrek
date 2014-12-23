#ifndef SNPESTIMATION_H
#define SNPESTIMATION_H

#include <vector>
#include "decomposition.h"
#include "genotypefilehandler.h"

class SnpEstimation
{
	public:
		SnpEstimation(GenotypeFileHandler *genotypeFileHandler,SnpIndex *snpIndex, std::vector<Snp*> *snpList);
		virtual ~SnpEstimation();
	protected:
	private:
		GenotypeFileHandler *m_genotypeFileHandler;
		SnpIndex *m_snpIndex;
		std::vector<Snp*> *m_snpList;
};

#endif // SNPESTIMATION_H
