#include "snpestimation.h"

SnpEstimation::SnpEstimation(const Command &commander){
    m_thread = commander.getNThread();
    m_ldCorrection = commander.ldCorrect();
    m_extreme = commander.getExtremeAdjust();
    m_prevalence = commander.getPrevalence();
    m_qt = commander.quantitative();
    m_bt = commander.binary();
}

SnpEstimation::~SnpEstimation()
{
    //dtor
}

void SnpEstimation::estimate(GenotypeFileHandler &genotypeFileHandler,const std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> &regionList){
    // This contains the genotypes from the reference panel, the main container for this function
    boost::ptr_deque<Genotype> genotype;
    // This contains the index of the SNPs included in the genotype
    // Important for everything (otherwise, how do we the beta and update the heritability?)
    std::deque<size_t> snpLoc;
    // This should be invoked when we came to the end of chromosome / region
    // Such that the buffed variance can be added back to the total variance
    // Useless when the sign of the statistics is not given
    bool finalizeBuff = false;
    fprintf(stderr, "Estimate SNP heritability\n");
    // Ignore the progress bar just yet
    // Only add it in when everything is completed
    /** Progress Bar related code here **/
    // The completed flag should be the return from GenotypeFileHandler. When all the SNPs
    // from the reference panel were used, then the completed flag should change to true
    bool completed = false;
    std::deque<size_t> boundary; //The boundary is used to indicate the blocks
    while(!completed){
        // Get the required genotypes
        genotypeFileHandler.getSnps(snpIndex, snpList, genotype, snpLoc, finalizeBuff, completed, boundary);

    }
}
