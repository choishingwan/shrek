#ifndef COMMAND_H
#define COMMAND_H

#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include "usefulTools.h"

class Command
{
    public:
        Command(int argc, char* argv[], bool &error);
        virtual ~Command();
        void printRunSummary(std::string regionMessage);
        void printBriefUsage();
        size_t Getthread() const;
        size_t GetminBlock() const;
        size_t GetmaxBlock() const;
        size_t GetsampleSize() const;
        size_t GetcaseSize() const;
        size_t GetcontrolSize() const;
        size_t GetcIndex() const;
        size_t GettIndex() const;
        size_t GetbpIndex() const;
        size_t GetchrIndex() const;
        size_t GetrsIndex() const;
        size_t GetsampleSizeIndex() const;
        double Getprevalence() const;
        double Getmaf() const;
        bool ldCorrect() const;
        bool validate() const;
        bool isPvalue() const;
        bool provideSampleSize() const;
        bool quantitative() const;
        bool caseControl() const;
        bool maxBlockSet() const;
        std::string GetoutputPrefix() const;
        std::string GetpValueFileName() const;
        std::string GetldFilePrefix() const;
        std::string GetregionList() const;
        std::string GetprogrammeName() const;
    protected:
    private:
        size_t m_thread;
        size_t m_minBlock;
        size_t m_maxBlock;
        size_t m_sampleSize;
        size_t m_caseSize;
        size_t m_controlSize;
        size_t m_cIndex;
        size_t m_tIndex;
        size_t m_bpIndex;
        size_t m_chrIndex;
        size_t m_rsIndex;
        size_t m_sampleSizeIndex;
        double m_prevalence;
        double m_maf;
        bool m_ldCorrection;
        bool m_validate;
        bool m_isPvalue;
        bool m_provideSampleSize;
        bool m_quantitative;
        bool m_caseControl;
        bool m_maxBlockSet;
        std::string m_outputPrefix;
        std::string m_pValueFileName;
        std::string m_ldFilePrefix;
        std::string m_regionList;
        std::string m_programmeName;
        void printUsage();
};

#endif // COMMAND_H
