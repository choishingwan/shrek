#include "snpestimation.h"

SnpEstimation::SnpEstimation(const Command &commander){
    m_thread = commander.getNThread();
    m_ldCorrection = commander.ldCorrect();
    m_extreme = commander.getExtremeAdjust();
    m_prevalence = commander.getPrevalence();
    m_qt = commander.quantitative();
    m_bt = commander.binary();
    m_blockSize = commander.getSizeOfBlock();
    m_output = commander.getOutputPrefix();
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
    bool starting = true;
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
            bool skip = false;
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
                skip = true; // don't bother in padding
            }
            else{} //everything as normal
            // now pad the remaining SNPs for the last block
            // then construct the LD and remove the perfectLD again.
            // however, this time, we are certain there will be no boundary change
            if(!skip){
                genotypeFileHandler.getSNP(snpList, genotype, snpLoc, finalizeBuff, completed,boundary);
                linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
            }
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
        if(finalizeBuff) decompose.run(linkage, snpLoc, boundary, snpList, finalizeBuff, !retainLastBlock, starting, regionList);
        else decompose.run(linkage, snpLoc, boundary, snpList, false, false, starting, regionList);


        if(finalizeBuff){
            starting = true;
            // We need to check if we are going to store the last block
            if(retainLastBlock){
                boost::ptr_list<Genotype>::iterator genoIter = genotype.begin();
                size_t nRemoveElements =std::distance(snpLoc.begin(), boundary[1]);
                std::advance(genoIter,nRemoveElements);
                genotype.erase(genotype.begin(), genoIter);
                snpLoc.erase(snpLoc.begin(), boundary[1]);
                boundary.pop_front();
                linkage.clear(nRemoveElements);
            }
            else{
                //clean everything
                snpLoc.clear();
                boundary.clear();
                genotype.clear();
                linkage.clear();
            }
        }
        else{
            starting = false;
            //Also clean everything except the last block
            boost::ptr_list<Genotype>::iterator genoIter = genotype.begin();
            size_t nRemoveElements =std::distance(snpLoc.begin(), boundary[1]);
            std::advance(genoIter,nRemoveElements);
            genotype.erase(genotype.begin(), genoIter);
            snpLoc.erase(snpLoc.begin(), boundary[1]);
            boundary.pop_front();
            linkage.clear(nRemoveElements);

        }
        fprintf(stderr, "Next round\n");
    }
    fprintf(stderr, "Estimated the SNP Heritability, now output the results\n");
}



void SnpEstimation::result(const boost::ptr_vector<Snp> &snpList, const boost::ptr_vector<Region> &regionList){

    double i2 = (m_bt)? usefulTools::dnorm(usefulTools::qnorm(m_prevalence))/(m_prevalence): 0;
    i2 = i2*i2;
    double adjust = 0.0;
    size_t count  = 0;
    double sampleSize = 0.0;
    std::vector<double> heritability;
    std::vector<double> variance;
    std::vector<double> effective;
    for(size_t i = 0; i < regionList.size(); ++i){
        // add stuff
        heritability.push_back(0.0);
        variance.push_back(regionList[i].getVariance());
        effective.push_back(0.0);
    }
    std::ofstream fullOutput;
    bool requireFullOut = false;
    if(!m_output.empty()) requireFullOut = true;
    if(requireFullOut){
        std::string fullOutName = m_output;
        fullOutName.append(".res");
        fullOutput.open(fullOutName.c_str());
        if(!fullOutput.is_open()){
            requireFullOut = false;
            fprintf(stderr, "Cannot access output file: %s\n", fullOutName.c_str());
            fprintf(stderr, "Will only provide basic output\n");
        }
        else{
            fullOutput << "Chr\tLoc\trsID\tZ\tH\tStatus" << std::endl;
        }
    }
    for(size_t i = 0; i < snpList.size(); ++i){
        if(requireFullOut){
            fullOutput << snpList[i].getChr() << "\t" << snpList[i].getLoc() << "\t" << snpList[i].getRs() << "\t" << snpList[i].getStat() << "\t" << snpList[i].getHeritability() << "\t" << snpList[i].getStatus() << std::endl;
        }
        for(size_t j = 0; j < regionList.size(); ++j){
            if(j ==0) std::cerr << snpList[i].flag(j) << std::endl;
            if(snpList[i].flag(j)) heritability[j]+=snpList[i].getHeritability();
            if(snpList[i].flag(j)) effective[j] += snpList[i].getEffective();
        }
        sampleSize +=(double)(snpList[i].getSampleSize());
        count++;
        if(m_bt){
            double portionCase = (double)(snpList[i].getNCase()) / (double)(snpList[i].getSampleSize());
            adjust+= ((1.0-m_prevalence)*(1.0-m_prevalence))/(i2*portionCase*(1-portionCase));
        }

    }

    if(requireFullOut) fullOutput.close();

    double adjustment = adjust/(double)count; // we use the average adjustment value here
    if(!m_bt) adjustment = m_extreme;
    double averageSampleSize = sampleSize/(double)count;
    std::cerr << "Adjustment is: " << adjustment << std::endl;
    if(!m_output.empty()) requireFullOut =true;
    std::ofstream sumOut;
    if(requireFullOut){
        std::string sumOutName = m_output;
        sumOutName.append(".sum");
        sumOut.open(sumOutName.c_str());
        if(!sumOut.is_open()){
            fprintf(stderr, "Cannot access output file %s\n", sumOutName.c_str());
            fprintf(stderr, "Will output to stdout instead\n");
            requireFullOut=false;
        }
    }
    if(!requireFullOut){
        // Only output to the stdout
        std::cout << "Region\tHeritability\tVariance" << std::endl;
        for(size_t i = 0; i < regionList.size(); ++i){
            std::cout << regionList[i].getName() << "\t" << adjustment*heritability[i] << "\t";

            if(variance[i]==0.0){ // We use effective number
                std::cout << adjustment*(2.0*(effective[i]*adjustment*adjustment+2.0*adjustment*effective[i]*averageSampleSize)/(averageSampleSize*averageSampleSize)) << std::endl;
            }
            else{
                std::cout << adjustment*adjustment*variance[i] << std::endl;
            }
        }
    }
    else{
        sumOut << "Region\tHeritability\tVariance" << std::endl;
        for(size_t i = 0; i < regionList.size(); ++i){
            sumOut << regionList[i].getName() << "\t" << adjustment*heritability[i] << "\t";

            if(variance[i]==0.0){ // We use effective number
                sumOut << adjustment*(2.0*(effective[i]*adjustment*adjustment+2.0*adjustment*effective[i]*averageSampleSize)/(averageSampleSize*averageSampleSize)) << std::endl;
            }
            else{
                sumOut << adjustment*adjustment*variance[i] << std::endl;
            }
        }
        sumOut.close();
    }
}
