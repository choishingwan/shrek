#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <cstdlib>
#include <string.h>
#include <iostream>
#include <deque>

class Genotype
{
	public:
		Genotype();
		virtual ~Genotype();
        double Getr(Genotype* snpB, bool correction);
        double GetrSq(Genotype* snpB, bool correction);
        static void SetsampleNum(size_t sampleNum);
        void Setmean(double mean);
        void SetstandardDeviation(double standardDeviation);
        void AddsampleGenotype(int genotype, size_t sampleIndex);
        static void clean(std::deque<Genotype*> &genotype, size_t remaining);
	protected:
	private:
        unsigned long long *m_genotypeA;
        unsigned long long *m_genotypeB;
        static size_t m_sampleNum;
        unsigned int m_bitSize;
        unsigned int m_requiredBit;
		double m_mean;
		double m_standardDeviation;
        unsigned int m_sum; //Might want to check if it will go out of bound
};

#endif // GENOTYPE_H
