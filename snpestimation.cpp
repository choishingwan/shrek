#include "snpestimation.h"

SnpEstimation::SnpEstimation(const Command &commander){
    m_thread = commander.getNThread();
    m_ldCorrection = commander.ldCorrect();
    m_extreme = commander.getExtremeAdjust();
    m_prevalence = commander.getPrevalence();
    m_qt = commander.quantitative();
    m_bt = commander.binary();
    m_blockSize = commander.getSizeOfBlock();
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

    Linkage linkage(m_thread, m_blockSize);
    Decomposition decompose(m_thread);
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
    size_t roundNumber = 0;
    while(!completed){
        // Get the required genotypes
        // The concept now is simpler
        // 1. Get all the required SNPs for the analysis using genotypeFileHandler
        //    - Remove all SNPs that doesn't pass the threshold or QC
        // 2. Perform the LD matrix construction
        // 3. For each window, remove the perfect LD by not including them in the decomposition
        //    - Remember, if they are in perfect LD, we can just move the row around and they will still be the same
        fprintf(stderr, "Get Block\n");
        bool retainLastBlock=false; // only used when the finalizeBuff is true, this indicate whether if the last block is coming from somewhere new
        while(boundary.size() < 4 && !completed && !finalizeBuff){
            // We will continue to get block
            // When getBlock return, the iterator are always pointing to the next valid SNP
            // However, it is possible for this SNP to not pass the MAF threshold or are in
            // perfect LD with SNPs in previous blocks. Therefore we need to make sure that
            // when we remove/filter this SNP, the next SNP is still within the acceptable
            // region
            genotypeFileHandler.getBlock(snpList, genotype, snpLoc, finalizeBuff, completed,boundary);
            // Now calculate the LD
            bool boundChange = false;
            linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
            // Three possibilities
            // 1. Boundary intact, then we can continue
            // 2. Boundary changed, but still in range. Read more SNPs and continue
            // 3. Boundary changed, out of range now, start decomposition
            // Here we need to check the boundaries
            if(boundChange && boundary.back() != snpLoc.end()){
                // need to check if the blocks are now ok, if not, change finalizeBuff to true
                // Also, need special
                size_t lastLoc = *(boundary.back());
                size_t prevLoc = *std::prev(boundary.back());
                if(snpList.at(lastLoc).getLoc()-snpList.at(prevLoc).getLoc() > m_blockSize){
                    // this is problematic here
                    retainLastBlock = true;
                    finalizeBuff = true;
                }
                else{}// everything as normal

            }
            else if(boundChange){ // this is the abnormal situation where the after removing the perfect LDs, we lost the whole block
                finalizeBuff = true;
                boundary.pop_back(); //The last one is just for checking
            }
            else{} //everything as normal
            // now pad the remaining SNPs for the last block
            // then construct the LD and remove the perfectLD again.
            // however, this time, we are certain there will be no boundary change
            genotypeFileHandler.getSNP(snpList, genotype, snpLoc, finalizeBuff, completed,boundary);
            linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
        }
        // Need to check if we need to merge the last block into the one in front
        // Only consider it when we reached the end of the current region
        if(finalizeBuff && boundary.size() > 3){
//            // Now check if we need to merge the last two blocks
            size_t lastSnpOfThirdBlock = snpList.at(*std::prev(boundary.back())).getLoc();
            size_t lastSnp = snpLoc.back(); //The last of snpLoc is the last of the last block
            lastSnp = snpList.at(lastSnp).getLoc();
            if(lastSnp-lastSnpOfThirdBlock <= m_blockSize){
                boundary.pop_back(); //Because the boundary changed, we should also update the LD matrix
                // Mainly, the top right and bottom left corner
                bool boundChange = false;
                // Can consider writing a more efficient algorithm for this, however, this should also be rare.
                linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
                //There will be no boundary change because we have already removed the only bound that might be changed by the perfect LD
            }
        }
        // When we reach here, we are ready for decomposition
        fprintf(stderr, "Decomposition\n");
        if(finalizeBuff) decompose.run(linkage, snpLoc, boundary, snpList, !retainLastBlock, roundNumber, regionList);
        else decompose.run(linkage, snpLoc, boundary, snpList, false, roundNumber, regionList);

        linkage.print();

        if(finalizeBuff) roundNumber=0;
        else roundNumber++;

        // A number of possibilities
        // 1. Normal situation (1 last block that don't need to be checked)
        // 2. Ending (decompose the whole thing)
        // 3. Ending (1 last blcok extra)
        break; //End things first

    }
}
