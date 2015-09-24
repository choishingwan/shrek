// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef LINKAGE_H
#define LINKAGE_H

class LinkageThread;
#include <deque>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <thread>
#include <mutex>
//#include <unsupported/Eigen/IterativeSolvers>
#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <algorithm>
#include <limits>
#include <map>
#include <complex>
#include <math.h>
#include "configure.h"
#include "genotype.h"
#include "linkagethread.h"
#include "processcode.h"
#include "snp.h"

/** \class Linkage
 *  \brief Responsible for the LD matrix
 *
 *  This class is responsible for anything related to the LD matrix. Most importantly,
 *  it is responsible for the construction of LD matrix and the core decomposition step
 *  of the matrix equation.
 *  The most complicated part of this class is the part to deal with perfect LD. We take
 *  additional windows so that we make sure that for any window that were processing,
 *  their perfect LD partners are taken into account. This is based on the assumption that
 *  there will be no true perfect LD between two Snps if they are 1 window size away from
 *  each others.
 */
class Linkage
{
	public:
		Linkage(size_t thread);
		virtual ~Linkage();
        void Initialize(boost::ptr_deque<Genotype> &genotype, const size_t &prevResiduals);
        void Construct(boost::ptr_deque<Genotype> &genotype, const size_t &genotypeIndex, const size_t& remainedLD, const boost::ptr_vector<Interval> &blockSize, bool correction, std::deque<size_t> &ldLoc);
	    void print();
	protected:
	private:
	    void buildLd(bool correction, size_t vStart, size_t vEnd, size_t hEnd, boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &ldLoc);
        Eigen::MatrixXd m_linkage;
        size_t m_thread;
        static std::mutex mtx;

};

#endif // LINKAGE_H
