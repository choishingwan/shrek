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
        /** return the number of thread use */
        size_t Getthread() const;
        /** return the minimal block size requirement */
        size_t GetminBlock() const;
        /** return the maximum block size restriction */
        size_t GetmaxBlock() const;
        /** return the sample size used */
        size_t GetsampleSize() const;
        /** return the number of cases */
        size_t GetcaseSize() const;
        /** return the number of control */
        size_t GetcontrolSize() const;
        /** return the column containing the test statistic/p-value for case control study */
        size_t GetcIndex() const;
        /** return the column containing the test statistic/p-value for quantitative trait study */
        size_t GettIndex() const;
        /** return the column containing the bp information */
        size_t GetbpIndex() const;
        /** return the column containing the chromosome information */
        size_t GetchrIndex() const;
        /** return the column containing the rs id information */
        size_t GetrsIndex() const;
        /** return the column containing the sample size information */
        size_t GetsampleSizeIndex() const;
        /** return the prevalence */
        double Getprevalence() const;
        /** return the maf threshold */
        double Getmaf() const;
        /** return the extreme adjustment parameter */
        double GetextremeAdjust() const;
        /** return whether if r-square correction is required */
        bool ldCorrect() const;
        /** return whether if snp validation is required */
        bool validate() const;
        /** return whether if p-value instead of test-statistic is given */
        bool isPvalue() const;
        /** return whether if sample size is defined */
        bool provideSampleSize() const;
        /** return whether if it is a quantitative trait study */
        bool quantitative() const;
        /** return whether if it is a case control study */
        bool caseControl() const;
        /** return whether if maximum block size is set */
        bool maxBlockSet() const;
        /** return whether if the p-value file contain a header */
        bool hasHeader() const;
        /** return output prefix */
        std::string GetoutputPrefix() const;
        /** return p-value file name */
        std::string GetpValueFileName() const;
        /** return genotype file prefix (for ld construction) */
        std::string GetldFilePrefix() const;
        /** return the list of region(s) information */
        std::string GetregionList() const;
        /** return the programme name */
        std::string GetprogrammeName() const;
        /** return the direction file name */
        std::string GetdirectionFile() const;
    protected:
    private:
        size_t m_thread; //!< Number of thread used
        size_t m_minBlock; //!< minimum block size required
        size_t m_maxBlock; //!< maximum block size allowed
        size_t m_sampleSize; //!< sample size
        size_t m_caseSize; //!< number of cases
        size_t m_controlSize; //!< number of control
        size_t m_cIndex; //!< the column containing the chi square information (or p-value), for case control, 1-based
        size_t m_tIndex; //!< the column containing the t-statistic information (or p-value), for quantitative traits, 1-based
        size_t m_bpIndex; //!< the column containing the location of the snp
        size_t m_chrIndex; //!< the column containing the chromosome information
        size_t m_rsIndex; //!< the column containing the rs-id
        size_t m_sampleSizeIndex; //!< the column containing the sample size information. Not use if sample size is specified
        double m_prevalence; //!< the prevalence
        double m_maf; //!< the maf filtering threshold. Snps with maf less then this threshold will be filtered out
        double m_extremeAdjust; //!< the extreme adjustment parameter
        bool m_ldCorrection;//!< whether if the r square / r should be adjusted (default true)
        bool m_validate; //!< if the snp information should be validated
        bool m_isPvalue; //!< if the input is p-value instead of test statistics
        bool m_provideSampleSize;//!< whether if sample size information is provided
        bool m_quantitative; //!< whether if it is quantitative trait. Mutually exclusive with m_caseControl
        bool m_caseControl; //!< whether if it is case control. Mutually exclusive with m_quantitative
        bool m_maxBlockSet; //!< whether if the maximum block size is set
        bool m_providedMaf; //!< Indicate whether if the maf threshold is provided
        bool m_providedPrevalence; //!< Indicate whether if the prevalence information is provided
        bool m_provideExtremeAdjustment; //!< Indicate whether if the prevalence information is provided
        bool m_hasHeader; //!< Indicate whether if the p-value file contains header
        std::string m_outputPrefix; //!< the output prefix
        std::string m_pValueFileName; //!< the p-value input file
        std::string m_directionFile; //!< the direction file. Contain information of the direction of effect
        std::string m_ldFilePrefix; //!< the file prefix of the genotype file for the ld construction
        std::string m_regionList; //!< a list of region informations
        std::string m_programmeName; //!< the programme name. Use for the help message only
        /** Function to print the usage information of the programme */
        void printUsage();
        /** Function to perform general checking to for all parameters */
        bool generalCheck();
        /** Function to check parameters specific for case control studies */
        bool caseControlCheck();
        /** Function to check parameters specific for quantitative studies */
        bool quantitativeCheck();
};

#endif // COMMAND_H
