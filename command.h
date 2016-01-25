// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.

#ifndef COMMAND_H
#define COMMAND_H

#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <map>
#include "usefulTools.h"

/**
 * \class Command
 * \brief Parameter handler.
 *
 *  This is the class for parameter handling. It will hold all the basic user inputs
 *  and pass them to subsequent functions. It will also perform basic sanity check of
 *  the parameter inputs and contain the help messages
 *
 */


class Command
{
    public:
        /** Default constructor */
        Command();
        /** Default destructor */
        virtual ~Command();
        void initialize(int argc, char* argv[]);
        void printRunSummary(std::string regionMessage);
        void printBriefUsage();

        inline size_t getCaseSize() const{ return m_caseSize; };
        inline size_t getControlSize() const { return m_controlSize;};
        inline size_t getChr() const { return m_chrIndex-1; };
        inline size_t getRs() const { return m_rsIndex-1;};
        inline size_t getBp() const { return m_bpIndex-1;};
        inline size_t getSampleIndex() const {return m_sampleSizeIndex-1;};
        inline size_t getSign() const { return m_signIndex-1; };
        inline size_t getRef() const { return m_ref-1; };
        inline size_t getAlt() const { return m_alt-1; };
        inline size_t getImputeInfo() const {return m_imputeInfo-1;}

        inline size_t getThread() const { return m_thread;};
        inline size_t getDistance() const {return m_distance;};
        inline size_t getSampleSize() const {return m_sampleSize;};
        inline double getPrevalence() const { return m_prevalence;};
        inline double getMaf() const { return m_maf;};
        inline double getExtreme() const {return m_extremeAdjust;};
        inline bool validate() const { return m_validate;};
        inline bool isPvalue() const {return m_isPvalue;};
        inline bool ldCorrect() const {return m_ldCorrection;};
        inline bool quantitative() const {return m_qt;};
        inline bool caseControl() const {return m_cc;};
        inline bool conRisk() const {return m_rqt;};
        inline bool diRisk() const { return m_rcc;};
        inline bool mafFilter() const {return m_providedMaf; };
        inline bool extremeAdjust() const {return m_provideExtremeAdjustment;};
        inline bool sampleSizeProvided() const{return m_provideSampleSize;};
        inline bool removeAmbig() const{return m_keepAmbiguous; }; //Default keeping ambiguous SNPs
        inline bool hasSign() const{return m_signGiven; };
        inline std::string getPvalueFileName() const{return m_pValueFileName; };
        inline std::string getLdFilePrefix() const{ return m_ldFilePrefix; };
        inline std::string getOutputPrefix() const {return m_outputPrefix; };
        inline std::string getRegion() const{ return m_regionList;};
        inline std::string getGenotype() const{ return m_genotypeFilePrefix; };
        inline size_t getStat() const {return m_stats; };

    protected:
    private:
        /**
                Meta information
        */
        double m_version=0.03;
        // We will use the predefined programme name (SHREK)

        //Other information
        //Basically avoid any default values for the index
        size_t m_caseSize=0;
        size_t m_controlSize=0;
        size_t m_chrIndex = 0;
        size_t m_rsIndex = 0;
        size_t m_bpIndex = 0;
        size_t m_sampleSizeIndex = 0;
        size_t m_imputeInfo=0;
        size_t m_signIndex = 0; // Index for the direction of effect
        size_t m_ref = 0;
        size_t m_alt = 0;
        size_t m_thread = 1;
        size_t m_distance = 1000000;
        size_t m_sampleSize = 0;
        double m_prevalence=1.0;
        double m_maf = -1.0;
        double m_extremeAdjust = 1.0;
        double m_imputeThreshold=0.9;
        double m_nullSign = 0;
        bool m_validate=true; // Default run the validation check
        bool m_isPvalue = false;
        bool m_ldCorrection=true; //Always use the LD correction as the default
        bool m_qt = false;
        bool m_cc = false;
        bool m_rqt = false;
        bool m_rcc = false;
        bool m_providedPrevalence = false;
        bool m_providedMaf = false;
        bool m_provideExtremeAdjustment = false;
        bool m_provideSampleSize = false;
        bool m_keepAmbiguous = false; //Default not keeping ambiguous SNPs
        bool m_signGiven = false; // Whether if the direction is given

        std::string m_pValueFileName="";
        std::string m_ldFilePrefix="";
        std::string m_outputPrefix="";
        std::string m_regionList="";
        std::string m_genotypeFilePrefix="";
        //std::vector<size_t> m_stats;
        size_t m_stats;

        void caseControlProcess(int argc, char* argv[]);
        void quantitativeProcess(int argc, char* argv[]);
        //void continuousRiskProcess(int argc, char* argv[]);
        //void dichotomusRiskProcess(int argc, char* argv[]);


        /** Function to print the usage information of the programme */
        void printUsage();
        void printCCUsage();
        void printQuantUsage();
        void printRiskCCUsage();
        void printRiskQtUsage();
        /** Function to perform general checking to for all parameters */
        bool generalCheck();

        //std::vector<size_t> processRange(std::string input);
        void getIndex(const std::vector<std::string> &index, const std::string &pvalueFileName, std::vector<size_t> &indexResult);
};

#endif // COMMAND_H
