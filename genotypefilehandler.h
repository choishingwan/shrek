#ifndef GENOTYPEFILEHANDLER_H
#define GENOTYPEFILEHANDLER_H


#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
#include <fstream>
#include <map>
#include <bitset> //For plink
#include <stdio.h>
#include <armadillo> //The matrix stuff
#include "usefulTools.h"
#include "snp.h"
#include "genotype.h"


class GenotypeFileHandler
{
    public:
        /** Default constructor */
        GenotypeFileHandler();
        /** Default destructor */
        virtual ~GenotypeFileHandler();
        void initialize(const Command &commander, const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList);
        void getSNP(const std::map<std::string, size_t> &snpIndex, boost::ptr_deque<Snp> &snpList,  boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, bool &finalizeBuff, bool &completed, std::deque<size_t> &boundary);
    protected:
    private:
        std::string m_genotypeFilePrefix="";
        size_t m_thread=1;
        size_t m_nRefSample = 0;
        double m_mafThreshold =0.0;

        // This is the file handler for the bed and bim file
        std::ifstream m_bedFile;
        std::ifstream m_bimFile;
        // Only used for log
        size_t m_nDuplicated= 0; //Check for duplication within the REFERENCE
        size_t m_nInvalid=0; // SNPs that have different information as the p-value file
        size_t m_nAmbig=0; // Number of SNPs that are ambiguous
        size_t m_nSnp=0;
        bool m_keepAmbiguous = false;
        bool m_include = false;
        bool openPlinkBinaryFile(const std::string s, std::ifstream & BIT);
        std::string m_prevChr = "";
        size_t m_prevLoc = 0;
        std::map<std::string, bool> m_duplicateCheck;
};

#endif // GENOTYPEFILEHANDLER_H
