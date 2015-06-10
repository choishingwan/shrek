// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <cstdlib>
#include <string.h>
#include <iostream>
#include <deque>

/** \class Genotype
 *  \brief Responsible for storing the genotype information and calculate the LD
 *
 *  Uses the Pearson correlation to calculate the linkage between two genotypes
 *  and can perform r/r-square correction based on Shieh et al (2010)
 */
class Genotype
{
	public:
	    /** Default constructor */
		Genotype();
		/** Default destructor */
		virtual ~Genotype();
        /** \brief Calculate the R
         *  \param [in] snpB, the other genotype
         *  \param [in] correction, whether if bias adjustment should be performed
         *  [out] The Pearson correlation between the two genotypes
         */
        double Getr(Genotype* snpB, bool correction);
        /** \brief Calculate the Rsq
         *  \param [in] snpB, the other genotype
         *  \param [in] correction, whether if bias adjustment should be performed
         *  [out] The square of the Pearson correlation between the two genotypes
         */
        double GetrSq(Genotype* snpB, bool correction);
        /** Setting the maximum number of samples used for LD construction */
        static void SetsampleNum(size_t sampleNum);
        /** Set the mean of the genotype */
        void Setmean(double mean);
        /** Set the standard deviation of the genotype */
        void SetstandardDeviation(double standardDeviation);
        /** Adding a new sample to this genotype */
        void AddsampleGenotype(int genotype, size_t sampleIndex);
        /** Cleaning all the pointers */
        static void clean(std::deque<Genotype*> &genotype, size_t remaining);
	protected:
	private:
        unsigned long long *m_genotypeA;
        unsigned long long *m_genotypeB;
        unsigned long long *m_missing;
        static size_t m_sampleNum;
        unsigned int m_bitSize;
        unsigned int m_requiredBit;
		double m_mean;
		double m_standardDeviation;
};

#endif // GENOTYPE_H
