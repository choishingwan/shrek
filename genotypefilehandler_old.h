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
#include <map>
#include <Eigen/Dense>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
#include "usefulTools.h"
#include "snp.h"
#include "genotype.h"
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
        GenotypeFileHandler();
        virtual ~GenotypeFileHandler();

        void initialize(const Command &commander, const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Interval> &blockInfo);
        void getSnps(boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::deque<size_t> &ldLoc, bool &chromosomeStart, bool &chromosomeEnd, size_t &prevResidual, boost::ptr_vector<Interval> &blockInfo);
        void getSnps(boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, std::deque<size_t> &ldLoc, bool &chromosomeStart, bool &chromosomeEnd, size_t &prevResidual, boost::ptr_vector<Interval> &blockInfo, const std::vector<int> &genoInclusion, Eigen::MatrixXd &betaInfo);
        inline size_t getSampleSize() const {return m_ldSampleSize; };
        inline size_t getSnpNumber() const {return m_finalSnpNumber; };
protected:
private:
        std::string m_genotypeFilePrefix="";
        std::string m_outPrefix="";
        std::ifstream m_bedFile;
        std::ifstream m_genoFile;
        std::vector<int> m_inclusion;
        size_t m_ldSampleSize=0;
        size_t m_inputSnp=0;
        size_t m_thread=1;
        size_t m_snpIter=0;
        size_t m_genoIter=0; //index for genoInclusion
        size_t m_blockLoc=0;
        size_t m_finalSnpNumber = 0;

        bool openPlinkBinaryFile(const std::string s, std::ifstream & BIT);
        void skipSnps(size_t const skipNum);
        void buildBlocks(std::string bimFileName, boost::ptr_vector<Interval> &blockInfo, size_t distance);

};

#endif // GENOTYPEFILEHANDLER_H