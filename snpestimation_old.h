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
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
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
		SnpEstimation();
		/** Default destructor */
		virtual ~SnpEstimation();
		/** Initialize the estimation */
		void Estimate(GenotypeFileHandler &genotypeFileHandler,const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, const Region& regionInfo, const Command &commander,boost::ptr_vector<Interval> &blockInfo);
		void Predict(GenotypeFileHandler &genotypeFileHandler,const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, const Region& regionInfo, const Command &commander,boost::ptr_vector<Interval> &blockInfo, const std::vector<int> &genoInclusion);
        void getResult(const Command &commander, const Region &region, const std::map<std::string, size_t> &snpIndex, const boost::ptr_vector<Snp> &snpList);
	protected:
	private:
        static inline void loadbar(size_t x, size_t n);
};

#endif // SNPESTIMATION_H
