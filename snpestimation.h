#ifndef SNPESTIMATION_H
#define SNPESTIMATION_H

#include <vector>
#include <deque>
#include <fstream>
#include <complex>
#include <iomanip>
#include "processcode.h"
#include "decomposition.h"
#include "genotypefilehandler.h"
#include "genotype.h"
#include "linkage.h"
#include "region.h"

#include "configure.h"//DEBUG

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
        double m_maf;
        double m_effective;
        bool m_correction;
        static inline void loadbar(size_t x, size_t n);
};

#endif // SNPESTIMATION_H
