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
//    boost::ptr_list<Genotype> genotype;
    boost::ptr_deque<Genotype> genotype;
    // This contains the index of the SNPs included in the genotype
    // Important for everything (otherwise, how do we the beta and update the heritability?)
//    std::list<size_t> snpLoc;
    std::deque<size_t> snpLoc;
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
    // This is use for indicating whether if the whole genome is read
    bool completed = false;
    std::deque<size_t> boundary;
    bool starting = true;
    size_t checking = 0; //DEBUG
    while(!completed){
        // Keep doing this until the whole genome is read
        bool retainLastBlock=false; // only used when the finalizeBuff is true, this indicate whether if the last block is coming from somewhere new
        while(boundary.size() < 4 && !completed && !finalizeBuff){
            // We want to get until the end or the end of the current region
//            fprintf(stderr, "Get Block now\n");
            // First, we get all the SNPs within the block
            genotypeFileHandler.getBlock(snpList, genotype, snpLoc, finalizeBuff, completed,boundary);
            size_t snpLocSizeBefore = snpLoc.size();

            bool boundChange = false; // This is to indicate whether if the bound SNP is in perfect LD
//            fprintf(stderr, "Construct\n");
            // Now construct the LD matrix and also remove the perfect LD
            linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
//            fprintf(stderr, "Got the linkage\n");
            if(snpLocSizeBefore == snpLoc.size()){
                // This mean there is no perfect LD, but we will still assert anyway just in case
                assert(!boundChange && "Error: Not possible");
            }
            else{ // Perfect LD was removed
                if(boundChange){
                    if(boundary.back()== snpLoc.size()){
                        // This mean the last block is in full LD with the previous block (unlikely though)
                        finalizeBuff=true;
                        boundary.pop_back();
                    }
                    else{
                        // We have changed the boundary, so we need to check if the new
                        // block is too far away
                        size_t lastSnpLocIndex = boundary.back();
                        size_t prevSnpLocIndex = boundary.back()-1;
                        if(snpList.at(snpLoc[lastSnpLocIndex]).getLoc()-
                           snpList.at(snpLoc[prevSnpLocIndex]).getLoc() > m_blockSize){
                            retainLastBlock = true;
                            finalizeBuff=true;
                        }
                        // We will patch it to make sure we are always starting
                        // with something correct for each round
                        genotypeFileHandler.getSNP(snpList, genotype, snpLoc, finalizeBuff, completed,boundary);
                        linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
                    }
                }
                else{
                    // Bound doesn't change, much easier to handle
                    // Just patch the blocks and rebuild the LD matrix
                    // We know that the bound will not change now
                    genotypeFileHandler.getSNP(snpList, genotype, snpLoc, finalizeBuff, completed,boundary);
                    linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
                }
            }
        }

        // Check if we need to merge the last two blocks
        if(finalizeBuff){
            if(retainLastBlock){
                // We do nothing here, the reason is that the implementation in linkage construct doesn't
                // take into consideration of the last block. So although the last block is considered as
                // too far away in this situation, the LD construction will not take special consideration
                // in that and might, in unlikely circumstances change the boundary SNP, therefore leading
                // to a problematic complication that we have not anticipated
            }
            else if(boundary.size() > 2){ // normal checking
                // We have enough blocks to perform merging (it is stupid to merge stuff if there is only two blocks,
                // but whatever)
                assert(boundary.back() > 0 && "The boundary is too small..." );
                size_t indexOfLastSnpSecondLastBlock = snpLoc.at(boundary.back()-1); // just in case
                size_t lastSnp = snpLoc.back();
                if(snpList.at(lastSnp).getLoc()-snpList.at(indexOfLastSnpSecondLastBlock).getLoc() <= m_blockSize){
                    boundary.pop_back();
                    // Because we change the boundary here, we want to reconstruct the LD such that it fits
                    bool boundaryChange = false;
                    linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundaryChange);
                    assert(!boundaryChange && "Error: Boundary should not change here");
                }
            }
        }

//        std::cerr << "Tidy up: " << retainLastBlock << " " << finalizeBuff << " " << starting << " " << boundary.size() << std::endl;
//        fprintf(stderr, "Decomposition\n");
        if(finalizeBuff) decompose.run(linkage, snpLoc, boundary, snpList, finalizeBuff, !retainLastBlock, starting, regionList);
        else decompose.run(linkage, snpLoc, boundary, snpList, false, false, starting, regionList);
//        fprintf(stderr, "Finish decompose\n");

        if(retainLastBlock){
            // Then we must remove everything except the last block
            // because finalizeBuff must be true here
            size_t update = boundary.back();
            genotype.erase(genotype.begin(), genotype.begin()+update);
            snpLoc.erase(snpLoc.begin(), snpLoc.begin()+update);
            size_t lastIndex =boundary.back();
            boundary.clear();
            boundary.push_back(lastIndex);
            linkage.clear(update);
            starting = true;
            finalizeBuff = false;
            retainLastBlock = false;
        }
        else{
            if(finalizeBuff){
                snpLoc.clear();
                boundary.clear();
                genotype.clear();
                linkage.clear();
                starting = true;
                finalizeBuff = false;
                retainLastBlock = false;
            }
            else{
                starting = false;
                size_t update=boundary[1];
                genotype.erase(genotype.begin(), genotype.begin()+update);
                snpLoc.erase(snpLoc.begin(), snpLoc.begin()+update);
                boundary.pop_front();
                linkage.clear(update);
            }
        }
        // We need to update the boundary to correctly adjust for the changes in the snpLoc
        size_t updateNum = boundary.front();
        for(size_t i = 0; i < boundary.size(); ++i) boundary[i]-= updateNum;
//        linkage.print();
    }
    fprintf(stderr, "Estimated the SNP Heritability, now proceed to output\n");
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
            fullOutput << snpList[i].getChr() << "\t" << snpList[i].getLoc() << "\t" << snpList[i].getRs() << "\t" << snpList[i].getStat() << "\t" << snpList[i].getHeritability() << "\t"<< snpList[i].getStatus() << std::endl;
        }
        for(size_t j = 0; j < regionList.size(); ++j){
//            if(j ==0) std::cerr << snpList[i].flag(j) << std::endl;
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
//    std::cerr << "Adjustment is: " << adjustment << std::endl;
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
                std::cout << adjustment*(2.0*(effective[i]*adjustment*adjustment+2.0*adjustment*heritability[i]*averageSampleSize)/(averageSampleSize*averageSampleSize)) << std::endl;
            }
            else{
                std::cout << adjustment*adjustment*variance[i] << std::endl;
            }
        }
    }
    else{
        sumOut << "Region\tHeritability\tVariance" << std::endl;
        std::cout << "Region\tHeritability\tVariance" << std::endl;
        for(size_t i = 0; i < regionList.size(); ++i){
            sumOut << regionList[i].getName() << "\t" << adjustment*heritability[i] << "\t";
            std::cout << regionList[i].getName() << "\t" << adjustment*heritability[i] << "\t";

            if(variance[i]==0.0){ // We use effective number
                sumOut << adjustment*(2.0*(effective[i]*adjustment*adjustment+2.0*adjustment*heritability[i]*averageSampleSize)/(averageSampleSize*averageSampleSize)) << std::endl;
                std::cout << adjustment*(2.0*(effective[i]*adjustment*adjustment+2.0*adjustment*heritability[i]*averageSampleSize)/(averageSampleSize*averageSampleSize)) << std::endl;
            }
            else{
                sumOut << adjustment*adjustment*variance[i] << std::endl;
                std::cout << adjustment*adjustment*variance[i] << std::endl;
            }
        }
        sumOut.close();
    }
}
