// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef SNP_H
#define SNP_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <map>
#include "usefulTools.h"
#include "region.h"
#include "command.h"

/**
 * \class Snp
 * \brief Store Snp information and perform some basic logistic of these information
 */
class Snp
{
public:
        /** Default constructor */
        Snp(std::string chr, std::string rs, size_t bp, double sampleSize, double original, std::string refAllele, std::string altAllele);
        /** \brief Batch constructor for multiple Snps
         *  \param [out] snpList vector containing the pointer to all the Snps
         *  \param [in] commander structure that contains all the required
         *   instruction from command line parameters
         *
         *  This function is used to generate all the Snp information. It will
         *  fill up the Snp containers with the Snp informations and prepare it
         *  for use in the next function
         *  NOTE: The first line of the p-value file is considered to be the
         *  header information and will be skipped.
         */
        static void generateSnpList(std::vector<Snp*> &snpList, const Command *commander);
        /** \brief Generate Snp index, also remove duplicates (Quantitative trait)
         *
         *   This function is used for quantitative trait only.
         *  \param [out] snpIndex, the container for the snp index
         *  \param [in]  snpList, the Snp container, containing all the snp information
         *  \param [in]  regionList, the region information
         *  \param [in]  isPvalue, indicate whether if the input is p-value only
         *  \param [in]  extremeRatio, the extreme adjustment ratio
         */
        static void generateSnpIndex(std::map<std::string, size_t> &snpIndex, std::vector<Snp*> &snpList, Region *regionList, bool isPvalue, double extremeRatio);
        /** \brief Generate Snp index, also remove duplicates (Quantitative trait)
         *
         *  This function is used for case control study only.
         *  \param [out] snpIndex, the container for the snp index
         *  \param [in]  snpList, the Snp container, containing all the snp information
         *  \param [in]  caseSize, the number of cases
         *  \param [in]  controlSize, the number of controls
         *  \param [in]  regionList, the region information
         *  \param [in]  isPvalue, indicate whether if the input is p-value only
         */
        static void generateSnpIndex(std::map<std::string, size_t> &snpIndex, std::vector<Snp*> &snpList, const size_t &caseSize, const size_t &controlSize, const double &prevalence, Region *regionList, bool isPvalue);
        static void generateSnpIndex(std::map<std::string, size_t> &snpIndex, std::vector<Snp*> &snpList);
        virtual ~Snp();
        std::string Getchr() const;
        std::string GetrsId() const;
        size_t Getbp() const;
        size_t GetregionSize() const;
        size_t GetperfectId () const;
        size_t GetblockInfo() const;
        double GetsampleSize() const;
        double Getoriginal() const;
        double Getbeta() const;
        double Getheritability() const;
        double GeteffectiveNumber() const;
        double GetsignedSqrtChiSq() const;
        double Getvariance() const;
        double GetsnpLDSC() const;
        bool Concordant(std::string chr, size_t bp, std::string rsId) const;
        bool GetFlag(size_t index) const;
        void Setheritability(double heritability);
        void SeteffectiveNumber(double effective);
        void setFlag(size_t index, bool value);
        void shareHeritability( Snp* i );
        void Setvariance(double i );
        void SetsnpLDSC(double i );
        void Setvariance(double const sigma, double const sigmaSquared, double const sigmaPowerThree, double const sigmaPowerFour );
        void SetadditionVariance(double i );
        void Setsign(int directionEffect);
        void SetblockInfo(size_t blockInfo);
        static void addDirection(std::map<std::string, size_t> &snpIndex, std::vector<Snp*> &snpList,std::string dirFile);
        /** \brief Function use to clean the pointers from the vector
         *  \param [in] snpList, the vector containing the snp pointers
         */
        static void cleanSnp(std::vector<Snp*> &snpList);
        static bool sortSnp (Snp* i, Snp* j);
        /** The maximum sample size */
        static size_t GetmaxSampleSize();
        static double Getadjustment();
        static void Setadjustment(const double prevalence, const size_t caseSize, const size_t controlSize);


protected:
private:
        std::string m_chr;
        std::string m_rs;
        std::string m_ref;
        std::string m_alt;
        size_t m_bp;
        size_t m_blockInfo;
        int m_sign;
        double m_sampleSize;
        double m_original;
        double m_effectiveNumber;
        double m_variance;
        double m_additionVariance;
        double m_sigmaPowerThree;
        double m_sigmaPowerFour;
        double m_snpLDSC;
        size_t m_perfectLdId;
        std::shared_ptr<double> m_beta; //Average of all Snps with perfect LD
        std::shared_ptr<double> m_sqrtChiSq; //Sqrt of the chiSquare, considered the sign.
        std::shared_ptr<double> m_heritability; //The master heritability
        std::vector<bool> m_regionFlag;



        /** This is use to form linked list between Snps that are in perfect LD */
        Snp* m_targetClass; //The master heritability


        void computeVarianceExplainedChi(const size_t &caseSize, const size_t &controlSize, const double &prevalence, bool isPvalue);
        void computeVarianceExplainedChi(bool isPvalue, double extremeRatio);

        static size_t m_perfectId;
        static double m_adjustment;
        static size_t m_maxSampleSize;

};

#endif // SNP_H
