#ifndef GENOTYPEFILEHANDLER_H
#define GENOTYPEFILEHANDLER_H

#include <boost/ptr_container/ptr_vector.hpp>
//#include <boost/ptr_container/ptr_list.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
#include <list>
#include <deque>
#include <stdexcept>
#include <fstream>
#include <map>
#include <bitset> //For plink
#include <stdio.h>
#include <assert.h>
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
        void getSNP(boost::ptr_vector<Snp> &snpList, boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, bool &windowEnd, bool &completed, std::vector<size_t> &boundary);
        void getBlock(boost::ptr_vector<Snp> &snpList, boost::ptr_deque<Genotype> &genotype, std::deque<size_t> &snpLoc, bool &windowEnd, bool &completed, std::vector<size_t> &boundary, bool addition);

    protected:
    private:
        std::string m_genotypeFilePrefix="";
        size_t m_nRefSample = 0;
        double m_mafThreshold =0.0;

        // This is the file handler for the bed and bim file
        std::ifstream m_bedFile;
        // Only used for log
        size_t m_nFilter = 0;
        size_t m_blockSize=0;
        bool m_keepAmbiguous = false;
        bool openPlinkBinaryFile(const std::string s, std::ifstream & BIT);


        // These are indicating the last USED SNP
        std::string m_prevChr = "";
        size_t m_prevLoc = 0, m_snpLoc=0;
        size_t m_nBytes=0;
        size_t m_nSnpSkipped=0;
        size_t m_nRequiredSnps=0;
        //Method number two, the method I used before, the inclusion vector
        std::vector<int> m_inclusion;
        std::vector<int> m_snpLineNumber;
        size_t m_snpIter=0; //This is for the iteration of inclusion
        size_t m_nSnp=0;
};

#endif // GENOTYPEFILEHANDLER_H
