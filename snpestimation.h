#ifndef SNPESTIMATION_H
#define SNPESTIMATION_H

#include <vector>
#include <deque>
#include <fstream>
#include "processcode.h"
#include "decomposition.h"
#include "genotypefilehandler.h"
#include "genotype.h"
#include "linkage.h"
#include "region.h"

class SnpEstimation
{
	public:
		SnpEstimation(GenotypeFileHandler *genotypeFileHandler,SnpIndex *snpIndex, std::vector<Snp*> *snpList, size_t thread, double maf, bool correction);
		void Estimate();
		void Getresult(std::string outputPrefix);
		virtual ~SnpEstimation();
	protected:
	private:
		GenotypeFileHandler *m_genotypeFileHandler;
		SnpIndex *m_snpIndex;
		std::vector<Snp*> *m_snpList;
		size_t m_thread;
		double m_check; //DEBUG
		double m_bug; //DEBUG
        double m_maf;
        double m_effective;
        bool m_correction;
};

#endif // SNPESTIMATION_H
