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
#include "usefulTools.h"

/**
 * \class Command
 * \brief Parameter handler.
 *
 *  This is the class for parameter handling. It will hold all the basic user inputs
 *  and pass them to subsequent functions. It will also perform basic sanity check of
 *  the parameter inputs
 *
 */


class Command
{
    public:
        /** Default constructor */
        Command();
        /** Default destructor */
        virtual ~Command();
        /** initialize the command handler. read from the command line and get the parameters */
        void initialize(int argc, char* argv[]);
        /** Print the run summary */
        void printRunSummary(std::string regionMessage);
        /** Print the brief usage information */
        void printBriefUsage();
        inline size_t getCaseSize() const{ return m_caseSize; };

    protected:
    private:
        /**
                Meta information
        */
        double m_version;
        std::string m_programmeName; //!< the programme name. Use for the help message only

        //Other information
        size_t m_caseSize;
        size_t m_controlSize;
        size_t m_chrIndex = 0;
        size_t m_rsIndex = 1;
        size_t m_bpIndex = 2;
        size_t m_sampleSizeIndex = 3;
        size_t m_dirIndex = 7;
        size_t m_ref = 0;
        size_t m_alt = 0;
        size_t m_thread = 1;
        size_t m_maxBlock = 0;
        size_t m_minBlock = 0;
        size_t m_distance = 3000000;
        size_t m_sampleSize = 0;
        double m_prevalence=1.0;
        double m_maf = -1.0;
        double m_extremeAdjust = 1.0;
        bool m_maxBlockSet=false;
        bool m_validate=false;
        bool m_isPvalue = false;
        bool m_ldCorrection=false;
        bool m_qt = false;
        bool m_cc = false;
        bool m_rqt = false;
        bool m_rcc = false;
        bool m_providedPrevalence = false;
        bool m_providedMaf = false;
        bool m_provideExtremeAdjustment = false;
        bool m_provideSampleSize = false;
        bool m_keep = true; //Default keeping ambiguous SNPs

        std::string m_pValueFileName="";
        std::string m_ldFilePrefix="";
        std::string m_outputPrefix="";
        std::string m_regionList="";
        std::string m_genotypeFilePrefix="";
        std::vector<size_t> m_stats;

        void caseControlProcess(int argc, char* argv[]);
        void quantitativeProcess(int argc, char* argv[]);
        void continuousRiskProcess(int argc, char* argv[]);
        void dichotomusRiskProcess(int argc, char* argv[]);


        /** Function to print the usage information of the programme */
        void printUsage();
        void printCCUsage();
        void printQuantUsage();
        void printRiskUsage();
        /** Function to perform general checking to for all parameters */
        bool generalCheck();
        void riskMode(int argc, char* argv[]);
        void quantMode(int argc, char* argv[]);
        void ccMode(int argc, char* argv[]);
        std::vector<size_t> processRange(std::string input);
};

#endif // COMMAND_H
