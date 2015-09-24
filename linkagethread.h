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
#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

/** \class LinkageThread
 *  \brief The multi-threading component of the Linkage class
 *
 *  This class is responsible to handle the LD construction using
 *  multi-threading.
 */
class LinkageThread
{
	public:
        LinkageThread(bool correction, size_t vStart, size_t vEnd, size_t hEnd,  boost::ptr_deque<Genotype> *genotype, std::deque<size_t> *ldLoc,Eigen::MatrixXd *ldMatrix);
		virtual ~LinkageThread();

        static void *buildLd(void *in);
	protected:
	private:
        bool m_correction;
        size_t m_vStart;
        size_t m_vEnd;
        size_t m_hEnd;
        boost::ptr_deque<Genotype> m_genotype;
        std::deque<size_t> *m_ldLoc;
        Eigen::MatrixXd *m_ldMatrix;
        void ldConstruction();
        static std::mutex mtx;
};

#endif // LINKAGETHREAD_H
