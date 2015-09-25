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
#include <thread>
#include <mutex>
#include <Eigen/Dense>
#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "snp.h"
#include "region.h"
#include "linkage.h"
#include "genotype.h"

class Decomposition
{
	public:
		Decomposition(size_t thread);
        void run(const Linkage &linkage, const size_t &genotypeIndex, const size_t& remainedLD, const boost::ptr_vector<Interval> &blockInfo, std::deque<size_t> &ldLoc, const std::deque<size_t> &snpLoc, boost::ptr_vector<Snp> &snpList, const bool &chromosomeStart);
        void solve(const std::vector<size_t> boundaries, const size_t index, const Linkage &linkage, const std::deque<size_t> snpLoc, boost::ptr_vector<Snp> &snpList, const bool &chromosomeStart, const Eigen::MatrixXd &betaEstimate);
	protected:
	private:
        size_t m_thread=1;
        static std::mutex mtx;

};

#endif // DECOMPOSITION_H
