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
    boost::ptr_list<Genotype> genotype;
    // This contains the index of the SNPs included in the genotype
    // Important for everything (otherwise, how do we the beta and update the heritability?)
    std::list<size_t> snpLoc;
    // This should be invoked when we came to the end of chromosome / region
    // Such that the buffed variance can be added back to the total variance
    // Useless when the sign of the statistics is not given
    bool finalizeBuff = false;

    Linkage linkage(m_thread);
    fprintf(stderr, "Estimate SNP heritability\n");
    // Ignore the progress bar just yet
    // Only add it in when everything is completed
    /** Progress Bar related code here **/
    // The completed flag should be the return from GenotypeFileHandler. When all the SNPs
    // from the reference panel were used, then the completed flag should change to true
    bool completed = false;
    // This thing is getting advance...
    std::deque<std::list<size_t>::iterator > boundary;
    //std::deque<size_t> boundary; //The boundary is used to indicate the blocks
    while(!completed){
        // Get the required genotypes
        // The concept now is simpler
        // 1. Get all the required SNPs for the analysis using genotypeFileHandler
        //    - Remove all SNPs that doesn't pass the threshold or QC
        // 2. Perform the LD matrix construction
        // 3. For each window, remove the perfect LD by not including them in the decomposition
        //    - Remember, if they are in perfect LD, we can just move the row around and they will still be the same
        fprintf(stderr, "Get Block\n");
        genotypeFileHandler.getSNP(snpIndex, snpList, genotype, snpLoc, finalizeBuff, completed, boundary);

        // Now start performing the linkage stuff
        fprintf(stderr, "Compute LD\n");
        linkage.construct(genotype, snpLoc, boundary, snpList,m_ldCorrection);
        fprintf(stderr, "Done\n");
        break; //End things first

    }
}
