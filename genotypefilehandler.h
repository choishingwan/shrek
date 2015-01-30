#ifndef GENOTYPEFILEHANDLER_H
#define GENOTYPEFILEHANDLER_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <deque>
#include <string>
#include <bitset>
#include <math.h>
#include <fstream>
#include "usefulTools.h"
#include "snp.h"
#include "snpindex.h"
#include "genotype.h"
#include "processcode.h"
#include "genotype.h"

class GenotypeFileHandler
{
	public:
		GenotypeFileHandler(std::string genotypeFilePrefix, SnpIndex *snpIndex, std::vector<Snp*> &snpList, bool validate, bool maxBlockSet, size_t maxBlock, size_t minBlock, size_t thread);
		virtual ~GenotypeFileHandler();
        ProcessCode getSnps(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc, std::vector<Snp*> &snpList, bool &chromosomeStart, bool &chromosomeEnd, double const maf, size_t &prevResidual, size_t &blockSize);
        size_t GetsampleSize() const;
	protected:
	private:
		std::string m_genotypeFilePrefix;
        SnpIndex *m_blockSizeTract;
        SnpIndex *m_chrCount;
        size_t m_ldSampleSize;
        size_t m_expectedNumberOfSnp;
        size_t m_snpIter;
        size_t m_inputSnp;
        size_t m_thread;
        size_t m_processed;
        std::ifstream m_bedFile;
        bool openPlinkBinaryFile(const std::string s, std::ifstream & BIT);
        std::deque<std::string> m_chrExists;
        std::vector<int> m_inclusion;
        std::vector<size_t> m_locTract;

};

#endif // GENOTYPEFILEHANDLER_H
