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
        Snp(std::string chr, std::string rs, size_t bp, double sampleSize, double original, double beta);
        static void generateSnpList(std::vector<Snp*> &snpList, const std::string &pvalueFile, const size_t index, const size_t sampleSize, const size_t rsIndex, const size_t bpIndex, const size_t chrIndex, const size_t sampleIndex, bool sampleSizeProvided);
        static void generateSnpIndex(SnpIndex *snpIndex, std::vector<Snp*> &snpList, std::vector<std::vector<Region*> > &regionList, bool isPvalue, double extremeRatio);
        static void generateSnpIndex(SnpIndex *snpIndex, std::vector<Snp*> &snpList, const size_t &caseSize, const size_t &controlSize, const double &prevalence, std::vector<std::vector<Region*> > &regionList, bool isPvalue);
        virtual ~Snp();
        std::string Getchr() const;
        std::string GetrsId() const;
        size_t Getbp() const;
        size_t GetregionSize() const;
        double GetsampleSize() const;
        double Getoriginal() const;
        double Getbeta() const;
        double Getheritability() const;
        double GetheritabilityChi() const;
        //double Geteffective() const;
        bool Concordant(std::string chr, size_t bp, std::string rsId) const;
        bool GetFlag(size_t index) const;
        void Setheritability(double heritability);
        void setFlag(size_t index, bool value);
        void shareHeritability( Snp* i );
        //void Seteffective(double i );
        void Setvariance(double i );
        static void cleanSnp(std::vector<Snp*> &snpList);
        static bool sortSnp (Snp* i, Snp* j);
        static size_t m_maxSampleSize;
protected:
private:
        std::string m_chr;
        std::string m_rs;
        size_t m_bp;
        double m_sampleSize;
        double m_original;
        double m_oriBeta;
        double m_effectiveNumber;
        double m_variance;
        std::shared_ptr<double> m_beta; //Average of all Snps with perfect LD
        std::shared_ptr<double> m_heritability; //The master heritability
        Snp* m_targetClass; //The master heritability
        std::vector<bool> m_regionFlag;
        void computeVarianceExplainedChi(const size_t &caseSize, const size_t &controlSize, const double &prevalence, bool isPvalue);
        void computeVarianceExplained(const size_t &caseSize, const size_t &controlSize, const double &prevalence, bool isPvalue);
        void computeVarianceExplainedChi(bool isPvalue, double extremeRatio);
        void computeVarianceExplained(bool isPvalue);


};

#endif // SNP_H
