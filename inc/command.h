#ifndef COMMAND_H
#define COMMAND_H

#include <string>
#include <assert.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <map> //This is to check for duplicated columns
#include <stdio.h>
#include <iostream>
#include "usefulTools.h"

class Command
{
    public:
        /** Default constructor */
        Command();
        /** Default destructor */
        virtual ~Command();
        //Parse the command line arguments and prepare for downstream analysis
        bool parseCommand(int argc, char* argv[]);
        std::string getPValueFileName() const {return m_pValueFileName;};
        std::string getReferenceFilePrefix() const {return m_referenceFilePrefix;};
        std::string getRegion() const { return m_regionList; };
        size_t getMaxIndex() const {return m_maxIndex; };
        // Easiest way is to return stuff
        size_t getChrIndex() const {return m_chrIndex;};
        size_t getRsIndex() const {return m_rsIndex;};
        size_t getLocIndex() const {return m_bpIndex;};
        size_t getRefIndex() const {return m_refIndex;};
        size_t getAltIndex() const {return m_altIndex;};
        size_t getSampleIndex() const {return m_nSampleIndex;};
        size_t getImputeInfoIndex() const {return m_imputeInfoIndex;};
        size_t getSumStatIndex() const {return m_summaryStatisticIndex;};
        size_t getSignIndex() const {return m_signIndex;};
        size_t getCaseIndex() const {return m_caseIndex;};
        size_t getControlIndex() const {return m_controlIndex;};
        //Mode related
        bool quantitative() const {return m_qt;};
        bool binary() const {return m_caseControl;};
        // Other informations
        double getSignNull() const {return m_signNull;};
        int getNCase() const {return m_nCase;};
        int getNControl() const {return m_nControl;};
        double getPrevalence() const {return m_prevalence;};
        int getNSample() const {return m_nSample;};
        double getExtremeAdjust() const {return m_extremeAdjust;};
        double getMAF() const {return m_mafThreshold;};
        double getINFO() const {return m_infoThreshold;};
        int getNThread() const {return m_nThread;};
        int getSizeOfBlock() const {return m_sizeOfBlock;};
        bool isP() const {return m_isPvalue;};
        bool keepAmbiguous() const {return m_keepAmbiguous;};
        bool ldCorrect() const {return m_ldCorrect;};
        std::string getOutputPrefix() const {return m_outputPrefix;};



    protected:
    private:
        //metaInfo
        double m_version = 0.3;
        //Required Files
        std::string m_pValueFileName="";
        std::string m_referenceFilePrefix="";
        // File Column Index
        // 0 is reserved for null (parameter not provided)
        size_t m_maxIndex = 0; //This is for the size check when reading the p-value file
        size_t m_chrIndex = 0;
        size_t m_rsIndex = 0;
        size_t m_bpIndex = 0;
        size_t m_refIndex = 0;
        size_t m_altIndex = 0;
        size_t m_nSampleIndex = 0;
        size_t m_imputeInfoIndex=0;
        size_t m_summaryStatisticIndex = 0;
        size_t m_signIndex = 0;
        double m_signNull = 0.0;
        // Case Control
        bool m_caseControl = false;
        size_t m_caseIndex = 0;
        size_t m_controlIndex = 0;
        int m_nCase = 0;
        int m_nControl = 0;
        double m_prevalence = -1.0;

        // Quantitative Trait
        bool m_qt = false;
        int m_nSample = 0;
        size_t m_sampleIndex = 0;
        double m_extremeAdjust = 1.0;
        // Filtering information
        double m_mafThreshold = 0.0; // 0.0 = not performing filtering (nothing to filter)
        double m_infoThreshold = 0.0; // again, 0 = nothing to filter

        //Misc Information
        int m_nThread=1; // Use int for better checking of input
        int m_sizeOfBlock = 1000000; //Use int just in case this is overflow (e.g. -1);
        bool m_isPvalue = false; //Whether if the input is summary statistics or p-value
        bool m_keepAmbiguous = false;
        //bool m_validate = true; // Check whether if the SNPs on the p-value file is the same on the reference// Discourage no validation as that is a recipe for error
        bool m_ldCorrect = false; // Whether if the LD R2 should be corrected
        std::string m_regionList="";

        //Optional Files
        std::string m_outputPrefix="";

        //Functions
        bool processCode(int argc, char *argv[]);
        void usage();
        void btUsage();
        void qtUsage();
        void generalOptions();
};

#endif // COMMAND_H
