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
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/IterativeSolvers>
#include <algorithm>
#include <limits>
#include <map>
#include <complex>
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
	    /** The default constructor */
		Linkage();
		//Linkage(size_t thread, std::vector<Snp*> *snpList, std::deque<size_t> *snpLoc);
		/** The default destructor */
		virtual ~Linkage();
        /** \brief Construct the LD matrix
         *  \param [in] genotype, the genotype information
         *  \param [in] prevResidual, indicates how many have been done
         *  \param [in] blockSize, restriction on size of LD block
         *  \param [in] correction, whether if we want to correct the bias in R and Rsq
         */
		void Construct(std::deque<Genotype*> &genotype, const size_t &prevResidual, const size_t &blockSize, bool correction);
		/** Initialize the LD matrix */
		ProcessCode Initialize(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize);
		/** Reinitialize the matrix */
		ProcessCode Reinitialize(size_t &genotypeSize);
		size_t rows() const;
		size_t cols() const;
		/** Get the ld matrix of R square */
		Eigen::MatrixXd block(size_t blockStart, size_t lengthOfBlock);
		/** Get the ld matrix of R */
		Eigen::MatrixXd blockSqrt(size_t blockStart, size_t lengthOfBlock);
		/** Solving the matrix equation using the linkage matrix */
		Eigen::VectorXd solveChi(size_t start, size_t length, Eigen::VectorXd const *const betaEstimate, Eigen::VectorXd const *const sqrtChiSq, Eigen::MatrixXd *variance,Eigen::MatrixXd *additionVariance, size_t sampleSize);
		/** Removing perfectLD by setting them to 0 so that they will be updated */
        size_t Remove();
        /** Update the matrix after removing the perfect LD */
        void Update(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc);
        void setSnpLoc(std::deque<size_t> *snpLoc);
        void setSnpList(std::vector<Snp* > *snpList);
        void setThread(size_t thread);
		void print();

		static size_t DEBUG;
	protected:
	private:
        Eigen::MatrixXd m_linkage;
        Eigen::MatrixXd m_linkageSqrt;
        size_t m_thread;
        std::vector<size_t> m_perfectLd; //Store the remove index of on matrix level
        std::vector<Snp*> *m_snpList;
        std::deque<size_t> *m_snpLoc;
        void triangularThread(const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype);
		void rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype);
};

#endif // LINKAGE_H
