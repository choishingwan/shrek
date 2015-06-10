// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
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

/**\class GenotypeFileHandler
 * \brief Responsible for reading the plink genotype file
 *
 * This class is responsible to read enough Snps for the subsequent
 * process e.g. LD construction.
 * Currently we only support plink file input. In the future, we would
 * like to also support vcf file input.
 */
class GenotypeFileHandler
{
public:
        /** Default constructor */
        GenotypeFileHandler(std::string genotypeFilePrefix, size_t thread);
        void initialize(SnpIndex *snpIndex, std::vector<Snp*> *snpList, bool validate, bool maxBlockSet, size_t maxBlock, size_t minBlock);
        virtual ~GenotypeFileHandler();
        ProcessCode getSnps(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc, std::vector<Snp*> *snpList, bool &chromosomeStart, bool &chromosomeEnd, double const maf, size_t &prevResidual, size_t &blockSize);
        ProcessCode getSnps(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc, std::vector<Snp*> *snpList, bool &chromosomeStart, bool &chromosomeEnd, double const maf, size_t &numSnp);
        size_t GetsampleSize() const;
        size_t GetestimateSnpTotal() const;
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
        size_t m_estimateTotal;
        size_t m_defaultDistance;
        std::ifstream m_bedFile;
        bool openPlinkBinaryFile(const std::string s, std::ifstream & BIT);
        std::deque<std::string> m_chrExists;
        std::vector<int> m_inclusion;
        std::vector<size_t> m_locTract;

};

#endif // GENOTYPEFILEHANDLER_H
