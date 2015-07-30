// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef LINKAGETHREAD_H
#define LINKAGETHREAD_H

#include <deque>
#include <mutex>
#include <limits>
#include "linkage.h"
#include "configure.h"
#include "genotype.h"
#include "snp.h"
#include <math.h>

/** \class LinkageThread
 *  \brief The multi-threading component of the Linkage class
 *
 *  This class is responsible to handle the LD construction using
 *  multi-threading.
 */
class LinkageThread
{
	public:
		LinkageThread(bool correction, const size_t blockEnd, Eigen::MatrixXd *ldMatrix, Eigen::MatrixXd *ldMatrixSqrt, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList);
		LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, Eigen::MatrixXd *ldMatrix, Eigen::MatrixXd *ldMatrixSqrt, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList);
		LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd , const size_t boundEnd, Eigen::MatrixXd *ldMatrix, Eigen::MatrixXd *ldMatrixSqrt, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList);
		virtual ~LinkageThread();

		void Addstart(size_t i);
        static void *triangularProcess(void *in);
        static void *rectangularProcess(void *in);
        static void *simpleProcess(void *in);
	protected:
	private:
        bool m_correction;
        size_t m_snpStart;
        size_t m_snpEnd;
        size_t m_boundStart;
        size_t m_boundEnd;
        Eigen::MatrixXd *m_ldMatrix;
        Eigen::MatrixXd *m_ldMatrixSqrt;
        std::deque<Genotype*> *m_genotype;
        std::deque<size_t> *m_snpLoc;
        std::vector<size_t> m_startLoc;
        std::vector<size_t> *m_perfectLd;
        std::vector<Snp*> *m_snpList;
		void triangularProcess();
        void rectangularProcess();
        void simpleProcess();
        static std::mutex mtx;
};

#endif // LINKAGETHREAD_H
