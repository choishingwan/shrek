#ifndef GENOTYPEFILEHANDLER_H
#define GENOTYPEFILEHANDLER_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <deque>
#include <string>
#include <bitset>
#include "usefulTools.h"
#include "snp.h"
#include "snpindex.h"

class GenotypeFileHandler
{
	public:
		GenotypeFileHandler(std::string genotypeFilePrefix, SnpIndex *snpIndex, std::vector<Snp*> &snpList);
		virtual ~GenotypeFileHandler();
	protected:
	private:
		std::string m_genotypeFilePrefix;
        SnpIndex *m_chrCount;
        size_t m_ldSampleSize;
        size_t m_expectedNumberOfSnp;
        bool openPlinkBinaryFile(const std::string s, std::ifstream & BIT);
        std::vector<std::string> m_chrExists;
        std::vector<int> m_inclusion;
};

#endif // GENOTYPEFILEHANDLER_H
