#ifndef SNP_H
#define SNP_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <memory>
#include "snpindex.h"
#include "usefulTools.h"
#include "region.h"

class Snp
{
	public:
		Snp(std::string chr, std::string rs, size_t bp, size_t sampleSize, double original, double beta);
		static void generateSnpList(std::vector<Snp*> &snpList, const std::string &pvalueFile, const size_t index, const size_t sampleSize, const size_t rsIndex, const size_t bpIndex, const size_t chrIndex, const size_t sampleIndex, bool sampleSizeProvided);
		static void generateSnpIndex(SnpIndex *snpIndex, std::vector<Snp*> &snpList, std::vector<std::vector<Region*> > &regionList, bool isPvalue);
		static void generateSnpIndex(SnpIndex *snpIndex, std::vector<Snp*> &snpList, const size_t &caseSize, const size_t &controlSize, const double &prevalence, std::vector<std::vector<Region*> > &regionList, bool isPvalue);

		virtual ~Snp();
        std::string Getchr() const;
        std::string GetrsId() const;
        size_t Getbp() const;
        size_t GetsampleSize() const;
        size_t GetregionSize() const;
        double Getoriginal() const;
        double Getbeta() const;
        double GetDecomposeBeta() const;
        double Getheritability() const;
        double Geteffective() const;
        bool Concordant(std::string chr, size_t bp, std::string rsId) const;
        bool GetFlag(size_t index) const;
        void Setheritability(double heritability);
        void setFlag(size_t index, bool value);
<<<<<<< HEAD
        void shareHeritability(Snp* i);
=======
        void shareHeritability( Snp* i );
        void Seteffective(double i );
>>>>>>> perfectLd
        static void cleanSnp(std::vector<Snp*> &snpList);
	protected:
	private:
        std::string m_chr;
		std::string m_rs;
		size_t m_bp;
		size_t m_sampleSize;
        double m_original;
        double m_oriBeta;
<<<<<<< HEAD
        /** Self note:
				Normally we should also use the weak_ptr to prevent circular referencing. However, due to our algorithm, we should be able to avoid
				circular referencing. Therefore it will be much easier if we just use one instance of the pointer
         */
        std::shared_ptr<double> m_beta; //The master beta
        std::shared_ptr<size_t> m_betaCount; //The number of beta
        std::shared_ptr<double> m_heritability; //The master heritability

=======
        double m_effectiveNumber;
        std::shared_ptr<double> m_beta; //Average of all Snps with perfect LD
        std::shared_ptr<double> m_heritability; //The master heritability
>>>>>>> perfectLd
        std::vector<bool> m_regionFlag;
        static bool sortSnp (Snp* i, Snp* j);
        void computeVarianceExplained(const size_t &caseSize, const size_t &controlSize, const double &prevalence, bool isPvalue);
        void computeVarianceExplained(bool isPvalue);


};

#endif // SNP_H
