#ifndef RISKPREDICTION_H
#define RISKPREDICTION_H
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "snp.h"
#include "command.h"
#include "usefulTools.h"
#include "linkage.h"
#include "decomposition.h"
#include "genotypefilehandler.h"

class RiskPrediction
{
    public:
        /** Default constructor */
        RiskPrediction(const Command *commander,std::vector<Snp*> *snpList);
        /** Default destructor */
        virtual ~RiskPrediction();
        void checkGenotype();
        void run();
        void result();
    protected:
    private:
        std::vector<Snp*> *m_snpList;
        std::vector<bool> m_flipCheck;
        std::vector<int> m_genoInclude;
        std::vector<std::string> m_sampleId;
        std::vector<double> m_samplePheno;
        size_t m_thread;
        size_t m_minBlock;
        size_t m_maxBlock;
        double m_maf;
        bool m_ldCorrection;
        bool m_keep;
        bool m_maxBlockSet;
        bool m_validate;
        std::string m_genotypeFilePrefix;
        std::string m_ldFilePrefix;
        std::string m_outPrefix;
        std::map<std::string, size_t> snpIndex;
        GenotypeFileHandler *targetGenotype;
};

#endif // RISKPREDICTION_H
