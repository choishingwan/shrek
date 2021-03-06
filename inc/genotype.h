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
#include <limits.h>
#include <cmath>

/** \class Genotype
 *  \brief Responsible for storing the genotype information and calculate the LD
 *
 *  Uses the Pearson correlation to calculate the linkage between two genotypes
 *  and can perform r/r-square correction based on Shieh et al (2010)
 */
//typedef unsigned long long mlong;
typedef unsigned int mlong;

class Genotype
{
	public:
		Genotype();
		virtual ~Genotype();
        void GetbothR(const Genotype &snpB, const bool correction, double &r, double &rSq) const;
        static void SetsampleNum(size_t sampleNum);
        void AddsampleGenotype(const int first, const int second, const size_t sampleIndex);
        void Setmean(double mean);
        void SetstandardDeviation(double standardDeviation);
        size_t getNonMiss()const{return m_nonMissSample;};
	protected:
	private:
        mlong *m_genotype;
//        mlong *m_genotypeB;
        mlong *m_missing;
        double m_mean=0.0;
        double m_standardDeviation=0.0;
        size_t m_nonMissSample=0;
        static size_t m_sampleNum;
//        const static mlong m1 = 0x5555555555555555LLU;
//        const static mlong m2 = 0x3333333333333333LLU;
//        const static mlong m4 = 0x0f0f0f0f0f0f0f0fLLU;
//        const static mlong m8 = 0x00ff00ff00ff00ffLLU;
		const static mlong FIVEMASK = 0x55555555;
		const static mlong THREEMASK = 0x33333333;
		const static mlong OFMASK = 0x0f0f0f0f;
		const static mlong ONEZEROMASK = 0x01010101;
		const static mlong AAAAMASK = 0xaaaaaaaa;
//        const static mlong FIVEMASK = 0x5555555555555555LLU; //Not so simple
//		const static mlong THREEMASK = 0x3333333333333333LLU;
//		const static mlong AAAAMASK = 0xaaaaaaaaaaaaaaaaLLU;
//		const static mlong OFMASK = 0x0f0f0f0f0f0f0f0fLLU;
//		const static mlong ONEZEROMASK = 0x0101010101010101LLU;
        unsigned int m_bitSize;
        unsigned int m_requiredBit;

};
#endif // GENOTYPE_H
