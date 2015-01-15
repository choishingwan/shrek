#ifndef SNPESTIMATION_H
#define SNPESTIMATION_H

#include <vector>
#include <deque>
#include "processcode.h"
#include "decomposition.h"
#include "genotypefilehandler.h"
#include "genotype.h"

class SnpEstimation
{
	public:
		SnpEstimation(GenotypeFileHandler *genotypeFileHandler,SnpIndex *snpIndex, std::vector<Snp*> *snpList, size_t blockSize, size_t distance, size_t thread, double maf);
		void performEstimation();
		virtual ~SnpEstimation();
	protected:
	private:
		GenotypeFileHandler *m_genotypeFileHandler;
		SnpIndex *m_snpIndex;
		std::vector<Snp*> *m_snpList;
		size_t m_blockSize;
		size_t m_distance;
		size_t m_thread;
        double m_maf;
};

#endif // SNPESTIMATION_H
