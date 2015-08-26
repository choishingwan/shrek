// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include <vector>
#include <deque>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <map>
#include "snp.h"
#include "region.h"
#include "linkage.h"
#include "genotype.h"
#include "processcode.h"
#include "decompositionthread.h"

/** \class Decomposition
 *	\brief Responsible to generate the required input to solve matrix decomposition
 *
 *	This class is responsible for performing the decomposition of the equation
 *	RH=f
 *	Note: The actual solving happened within the linkage class
 */
class Decomposition
{
	public:
		/**Default constructor */
		Decomposition( std::vector<Snp*> *snpList, Linkage *linkageMatrix, size_t thread, Region *regionInfo);
		Decomposition( std::vector<Snp*> *snpList, Linkage *linkageMatrix, size_t thread);
		/**Default destructor */
		virtual ~Decomposition();
        /** The decomposition processor */
		ProcessCode Decompose(const size_t &blockSize, std::deque<size_t> &snpLoc, std::deque<Genotype*> &genotype, bool chromosomeStart, bool chromosomeEnd);
	protected:
	private:
        std::vector<Snp*> *m_snpList;
        Linkage *m_linkage;
		size_t m_thread;
		Region *m_regionInfo;

};

#endif // DECOMPOSITION_H
