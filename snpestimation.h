// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef SNPESTIMATION_H
#define SNPESTIMATION_H

#include <vector>
#include <deque>
#include <fstream>
#include <complex>
#include <iomanip>
#include <map>
#include "processcode.h"
#include "decomposition.h"
#include "genotypefilehandler.h"
#include "genotype.h"
#include "linkage.h"
#include "region.h"

#include "configure.h"//DEBUG

/** \class SnpEstimation
 *  \brief the driver class for Snp Heritability estimation
 *
 *  This class is responsible for calling and coordinate other functions
 *  that are required for the Snp calling.
 *
 *  In the future, if Risk prediction is also included, we can also then
 *  have another driver class without modify or use the SnpEstimation class
 */
class SnpEstimation
{
	public:
	    /** Default constructor */
		SnpEstimation(GenotypeFileHandler *genotypeFileHandler,std::map<std::string, size_t> *snpIndex, std::vector<Snp*> *snpList, size_t thread, double maf, bool correction, Region* regionInfo);
		/** Default destructor */
		virtual ~SnpEstimation();
		/** Initialize the estimation */
		void Estimate();
		/** Output the result of the estimation */
		void Getresult(std::string outputPrefix);
	protected:
	private:
		GenotypeFileHandler *m_genotypeFileHandler;
		std::map<std::string, size_t> *m_snpIndex;
		std::vector<Snp*> *m_snpList;
		size_t m_thread;
        double m_maf;
        double m_effective;
        bool m_correction;
        Region *m_regionInfo;
        static inline void loadbar(size_t x, size_t n);
};

#endif // SNPESTIMATION_H
