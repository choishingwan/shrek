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
#include "usefulTools.h"
#include "snp.h"
#include "snpindex.h"
#include "genotype.h"

class GenotypeFileHandler
{
	public:
		GenotypeFileHandler(std::string genotypeFilePrefix, SnpIndex *snpIndex, std::vector<Snp*> &snpList, bool validate);
		virtual ~GenotypeFileHandler();
        ProcessCode getSnps(std::deque<Genotype*> &genotype, size_t distance, size_t anchor, std::string anchorChr, std::deque<size_t> snpLoc, std::vector<Snp*> &snpList);
        size_t GetSampleSize() const;
	protected:
	private:
		std::string m_genotypeFilePrefix;
        SnpIndex *m_chrCount;
        size_t m_ldSampleSize;
        size_t m_expectedNumberOfSnp;
        size_t m_snpIter;
        size_t m_inputSnp;
        bool openPlinkBinaryFile(const std::string s, std::ifstream & BIT);
        std::vector<std::string> m_chrExists;
        std::vector<int> m_inclusion;
};

#endif // GENOTYPEFILEHANDLER_H
