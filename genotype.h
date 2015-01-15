#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <cstdlib>
#include <string.h>
#include <iostream>

class Genotype
{
	public:
		Genotype();
		virtual ~Genotype();
        double GetR(Genotype* snpB, bool correction);
        static void SetSampleNum(size_t sampleNum);
        void SetMean(double mean);
        void SetStandardDeviation(double standardDeviation);
        void AddSampleGenotype(int genotype, size_t sampleIndex);
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
