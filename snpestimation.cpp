#include "snpestimation.h"

SnpEstimation::SnpEstimation(const Command &commander){
    m_thread = commander.getNThread();
    m_ldCorrection = commander.ldCorrect();
    m_extreme = commander.getExtremeAdjust();
    m_prevalence = commander.getPrevalence();
    m_qt = commander.quantitative();
    m_bt = commander.binary();
    m_blockSize = commander.getSizeOfBlock(); // A better word might actually be the range of block?
    m_output = commander.getOutputPrefix();
}

SnpEstimation::~SnpEstimation()
{
    //dtor
}

void SnpEstimation::estimate(GenotypeFileHandler &genotypeFileHandler, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> &regionList){
    // This contains the genotypes from the reference panel, the main container for this function
    boost::ptr_deque<Genotype> genotype;
    // This contains the index of the SNPs included in the genotype
    std::deque<size_t> snpLoc;
    // This should be invoked when we came to the end of chromosome / region
    bool windowEnd = false;
    // Initialize the linkage and decompose
    Linkage linkage(m_thread);
    Decomposition decompose(m_thread);
    fprintf(stderr, "Estimate SNP heritability\n\n");

    size_t doneItems = 0;
    size_t totalSnp = snpList.size();
    std::string chr = "";

    // This is use for indicating whether if the whole genome is read
    bool completed = false;
    std::vector<size_t> boundary;
    bool starting = true;
    size_t checking = 0; //DEBUG
    while(!completed){
        // Keep doing this until the whole genome is read
        progress(doneItems, totalSnp, chr);
        // only used when the finalizeBuff is true, this indicate whether if the last block is coming from somewhere new
        bool retainLastBlock=false;
        while(boundary.size() < 5 && !completed && !windowEnd){
            genotypeFileHandler.getBlock(snpList, genotype, snpLoc, windowEnd, completed,boundary, false);
            bool boundChange = false;
            linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
            if(boundChange && boundary.back()==snpLoc.size()){
                    windowEnd=true;
                    boundary.pop_back();
            }
            else if(boundChange&&snpList.at(snpLoc[boundary.back()]).getLoc()- snpList.at(snpLoc[boundary.back()-1]).getLoc() > m_blockSize){
                    retainLastBlock = true;
                    windowEnd=true;
            }
            genotypeFileHandler.getBlock(snpList, genotype, snpLoc, windowEnd, completed,boundary, true);
            linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
        }
        bool merged = false;
        if(windowEnd && !retainLastBlock && boundary.size() > 2){ //Check whether if we need to merge the two blocks
            size_t indexOfLastSnpOfSecondLastBlock = snpLoc.at(boundary.back()-1); // just in case
            size_t lastSnp = snpLoc.back();
            if(snpList.at(lastSnp).getLoc()-snpList.at(indexOfLastSnpOfSecondLastBlock).getLoc() <= m_blockSize){
                boundary.pop_back();
                bool boundChange=false;
                linkage.construct(genotype, snpLoc, boundary, snpList, m_ldCorrection, boundChange);
                merged=true;
            }
        }

        if(retainLastBlock && !windowEnd) throw std::runtime_error("Impossible combination of windowEnd and retain last block!");
        decompose.run(linkage, snpLoc, boundary, snpList, windowEnd, !retainLastBlock, starting, regionList,merged);
        doneItems= snpLoc.at(boundary.back());
        chr = snpList[snpLoc[boundary.back()]].getChr();

        if(retainLastBlock){
            // Then we must remove everything except the last block
            // because finalizeBuff must be true here
            size_t update = boundary.back();
            genotype.erase(genotype.begin(), genotype.begin()+update);
            snpLoc.erase(snpLoc.begin(), snpLoc.begin()+update);
            boundary.clear();
            boundary.push_back(0);
            linkage.clear(update);
            starting = true;
            windowEnd = false;
            retainLastBlock = false;
        }
        else if(windowEnd){
            snpLoc.clear();
            boundary.clear();
            genotype.clear();
            linkage.clear();
            starting = true;
            windowEnd = false;
            retainLastBlock = false;
        }
        else{
            starting = false;
            size_t update=boundary[1];
            genotype.erase(genotype.begin(), genotype.begin()+update);
            snpLoc.erase(snpLoc.begin(), snpLoc.begin()+update);
            linkage.clear(update);
            for(size_t i = 0; i < boundary.size()-1; ++i) boundary[i] = boundary[i+1]-update;
            boundary.pop_back();
        }
    }
    progress(totalSnp, totalSnp, "");
    fprintf(stderr, "\n\nEstimated the SNP Heritability, now proceed to output\n");
}



void SnpEstimation::result(const boost::ptr_vector<Snp> &snpList, const boost::ptr_vector<Region> &regionList){
    bool sign = !(snpList.front().getSign()==0);
    double i2 = (m_bt)? usefulTools::dnorm(usefulTools::qnorm(m_prevalence))/(m_prevalence): 0;
    i2 = i2*i2;
    double adjust = 0.0;
    size_t countNum  = 0;
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
            fullOutput << "Chr\tLoc\trsID\tZ\tF\tH\tStatus" << std::endl;
        }
    }
    for(size_t i = 0; i < snpList.size(); ++i){
        if(requireFullOut){
            fullOutput << snpList[i].getChr() << "\t" << snpList[i].getLoc() << "\t" << snpList[i].getRs() << "\t" << snpList[i].getTStat() << "\t"<< snpList[i].getFStat() << "\t" << snpList[i].getHeritability() << "\t" << snpList[i].getStatus() << std::endl;
        }
        for(size_t j = 0; j < regionList.size(); ++j){
//            if(j ==0) std::cerr << snpList[i].flag(j) << std::endl;
            if(snpList[i].flag(j)) heritability[j]+=snpList[i].getHeritability();
            if(snpList[i].flag(j)) effective[j] += snpList[i].getEffective();
        }
        sampleSize +=(double)(snpList[i].getSampleSize());
        countNum++;
        if(m_bt){
            double portionCase = (double)(snpList[i].getNCase()) / (double)(snpList[i].getSampleSize());
            adjust+= ((1.0-m_prevalence)*(1.0-m_prevalence))/(i2*portionCase*(1-portionCase));
        }

    }

    if(requireFullOut) fullOutput.close();

    double adjustment = adjust/(double)countNum; // we use the average adjustment value here
    if(!m_bt) adjustment = m_extreme;
    double averageSampleSize = sampleSize/(double)countNum;
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
    // We separate it, because it is possible for the previous if to change the requireFullOut flag
    if(!requireFullOut){
        // Only output to the stdout
        std::cout << "Region\tHeritability\tVariance" << std::endl;
        for(size_t i = 0; i < regionList.size(); ++i){
            std::cout << regionList[i].getName() << "\t" << adjustment*heritability[i] << "\t";
            if(!sign) std::cout << adjustment*(2.0*(effective[i]*adjustment*adjustment+2.0*adjustment*heritability[i]*averageSampleSize)/(averageSampleSize*averageSampleSize)) << std::endl;
            else std::cout << adjustment*adjustment*variance[i] << std::endl;
        }
    }
    else{
        sumOut << "Region\tHeritability\tVariance" << std::endl;
        std::cout << "Region\tHeritability\tVariance" << std::endl;
        for(size_t i = 0; i < regionList.size(); ++i){
            sumOut << regionList[i].getName() << "\t" << adjustment*heritability[i] << "\t";
            std::cout << regionList[i].getName() << "\t" << heritability[i] << "\t";
            if(!sign){ // We use effective number
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

void SnpEstimation::progress(size_t doneItems, size_t totalSnp, std::string chr){
    size_t barWidth = 50;
    double progress =(double)(doneItems)/(double)totalSnp;
    fprintf(stderr,"\r[");
    size_t pos=barWidth*progress;
    for(size_t i = 0; i < barWidth;++i){
        if(i < pos) fprintf(stderr,"=");
        else if(i==pos) fprintf(stderr,">");
        else fprintf(stderr," ");
    }
    fprintf(stderr, "]%u%% %s ", int(progress*100.0), chr.c_str() );
    fflush(stderr);
}
