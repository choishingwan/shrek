// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef DECOMPOSITIONTHREAD_H
#define DECOMPOSITIONTHREAD_H

#include <mutex>
#include <Eigen/Dense>
#include <deque>
#include <vector>
#include <iostream>
#include "linkage.h"
#include "snp.h"
#include "region.h"

/** \class DecompositionThread
 *  \brief Multi-threaded version of the decomposition.
 *
 *  This class is responsible to call the solve from linkage
 *  class to get the final result of the Snp heritability
 *  estimation
 */

class DecompositionThread
{
	public:

	    static Eigen::MatrixXd checking; //DEBUG

		/** Default constructor */
		DecompositionThread(size_t start, size_t length, Eigen::VectorXd const * const betaEstimate, Eigen::VectorXd const * const chiSq, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, Region *regionInfo);
		DecompositionThread(size_t start, size_t length, Eigen::MatrixXd const * const sampleMatrix,Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, std::vector<double> *samplePheno, bool chrStart);
		/** Default destructor */
		virtual ~DecompositionThread();
        /** Function called by the threading algorithm */
		static void *ThreadProcesser(void *in);
		static void *SampleProcesser(void *in);
	protected:
	private:
		size_t m_start;
		size_t m_length;
		size_t m_sampleSize;
		Eigen::VectorXd const * const m_betaEstimate;
		Eigen::VectorXd const * const m_sqrtChiSq;
		Eigen::MatrixXd const * const m_sampleMatrix;
		Linkage *m_linkage;
		std::deque<size_t> const *m_snpLoc;
        std::vector<Snp*> *m_snpList;
        std::vector<double> *m_samplePheno;
        bool m_chrStart;
        Region *m_regionInfo;
        static std::mutex decomposeMtx;

		/** Function to actually handle the solving */
		void solve();
		void sampleSolve();
        void fullProcess(Eigen::VectorXd const * const variance, Eigen::VectorXd const *const result, Eigen::VectorXd const *const effectiveReturnResult);
        void chromosomeStartProcess(Eigen::VectorXd const * const variance, Eigen::VectorXd const *const result, Eigen::VectorXd const *const effectiveReturnResult);
        void normalProcess(Eigen::VectorXd const * const variance, Eigen::VectorXd const *const result, Eigen::VectorXd const *const effectiveReturnResult);
        void endBlockProcess(Eigen::VectorXd const * const variance, Eigen::VectorXd const *const result, Eigen::VectorXd const *const effectiveReturnResult);

        void fullProcess(Eigen::MatrixXd const *const result, std::vector<double> *m_samplePheno);
        void chromosomeStartProcess(Eigen::MatrixXd const *const result, std::vector<double> *m_samplePheno);
        void normalProcess(Eigen::MatrixXd const *const result, std::vector<double> *m_samplePheno);
        void endBlockProcess(Eigen::MatrixXd const *const result, std::vector<double> *m_samplePheno);
};

#endif // DECOMPOSITIONTHREAD_H
