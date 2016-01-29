#include "snp.h"

Snp::Snp(std::string chr, std::string rsId, size_t loc, size_t nSample, size_t nCase, size_t nControl, std::string refAllele, std::string altAllele, double statistic, double info, int sign):m_chr(chr),m_rsId(rsId), m_loc(loc),m_nSample(nSample),m_nCase(nCase),m_nControl(nControl),m_ref(refAllele),m_alt(altAllele),m_statistic(statistic),m_infoScore(info),m_sign(sign){}

void Snp::setFlag(const size_t i, bool flag){ m_regionFlag.at(i) = flag;}
void Snp::setStatus(char status){ m_status=status; }
void Snp::flip(){
    std::string temp = m_ref;
    m_ref =m_alt;
    m_alt = temp;
    m_sign *= -1;
}

bool Snp::concordant(const std::string chr, const size_t loc, const std::string rs, std::string refAllele, std::string altAllele, bool &ambig){
    std::transform(refAllele.begin(), refAllele.end(), refAllele.begin(), toupper);
    std::transform(altAllele.begin(), altAllele.end(), altAllele.begin(), toupper);
    std::transform(m_ref.begin(), m_ref.end(), m_ref.begin(), toupper);
    std::transform(m_alt.begin(), m_alt.end(), m_alt.begin(), toupper);
    /** An important point here **/
    /** If we don't need to flip, we will never consider the SNP as ambiguous */
    //Transform everything to capital to avoid problem
    if(m_chr.compare(chr)==0 && m_loc == loc && m_rsId.compare(rs)==0){
        // Check if the reference and alternatives are provided
        // If either of them are empty, then we will just take it as if it is not ambiguous
        if(m_ref.empty() || m_alt.empty())  return true;
        else{
            if(m_ref.compare(refAllele)==0 && m_alt.compare(altAllele)==0) return true;
            else if(m_ref.compare(altAllele)==0 && m_alt.compare(refAllele)==0){
                //Flipping works, so we should flip
                flip();
                // In the following cases, we consider it to be ambiguous
                if( (m_ref.compare("A")==0 && m_alt.compare("T")==0) ||
                        (m_ref.compare("T")==0 && m_alt.compare("A")==0) ||
                        (m_ref.compare("G")==0 && m_alt.compare("C")==0) ||
                        (m_ref.compare("C")==0 && m_alt.compare("G")==0))
                    ambig=true;

                return true;
            }
            return false;
        }
        return true;
    }
    return false;
}
void Snp::computeSummaryStatistc(){
    double beta = 0.0;
    if(m_statistic < 1.0){ //Direct transform p-value of 1 to summary statistics of 0
        errno = 0;
        log(fabs(1.0-m_statistic/2.0)); //We know the reason of calculation problem is the log, so we check if the log will cause any error
        if(errno != 0){ //It means we cannot get the result from it
            errno = 0;
            log(fabs(m_statistic/2.0));
            assert(errno!=0&&"I am uncertain why this error occurs");
            if(errno == 0) beta = fabs(usefulTools::qnorm((m_statistic+0.0)/2.0));
        }
        else beta = fabs(usefulTools::qnorm(1.0-((m_statistic+0.0)/2.0)));
    }
    if(m_nSample != 0) m_statistic = beta*beta; //Because we want to keep the statistic as the summary stat instead of p-value
    else m_statistic = beta; //Because we want to keep the statistic as the summary stat instead of p-value
}

Snp::~Snp()
{
    //dtor
}
void Snp::generateSnpIndex(std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> const &regionList){
    std::vector<size_t> regionIncrementationIndex(regionList.size(), 0); // The current index of all regions are set to 0
    for(size_t i = 0; i < snpList.size(); ++i){
        snpIndex[snpList[i].getRs()] = i;
        // Now perform the flag setting
        snpList[i].m_regionFlag.push_back(false); //Set the base region flag to false;
        if(regionList.size() != 0){
            std::vector<bool> padding(regionList.size(), false);
            snpList[i].m_regionFlag.insert(snpList[i].m_regionFlag.end(), padding.begin(), padding.end()); //Initialize all the flags
        }
        for(size_t j = 1; j < regionList.size(); ++j){
            // Doesn't have to bother with the base region
            for(size_t k = regionIncrementationIndex.at(j); k < regionList.at(j).getIntervalSize(); ++k){
                if( regionList.at(j).getChr(k).compare(snpList[i].getChr())==0 &&
                    regionList.at(j).getStart(k) <= snpList[i].getLoc() &&
                    regionList.at(j).getEnd(k) >= snpList[i].getLoc()){
                    // There is one important assumption here:
                    // THE BED INTERVALS ARE NONOVERLAP
                    regionIncrementationIndex.at(j) = k;
                    snpList[i].setFlag(j, true);
                    break;
                }
            }
        }
    }
}



void Snp::generateSnpList(boost::ptr_vector<Snp> &snpList, const Command &commander){
    std::ifstream pValue;
    pValue.open(commander.getPValueFileName().c_str());
    if(!pValue.is_open()){
        throw std::runtime_error("Cannot read the p-value file");
    }
    std::string line;
    //Assume the p-value file has a header
    std::getline(pValue, line);
    bool qt = commander.quantitative();
    bool bt = commander.binary();
    bool isP = commander.isP();
    //If the index is 0, it means null
    size_t chrIndex = commander.getChrIndex();
    size_t rsIndex = commander.getRsIndex();
    size_t locIndex = commander.getLocIndex();
    size_t refIndex = commander.getRefIndex();
    size_t altIndex = commander.getAltIndex();
    size_t sampleIndex = commander.getSampleIndex();
    size_t imputeIndex = commander.getImputeInfoIndex();
    size_t sumStatIndex = commander.getSumStatIndex();
    size_t signIndex = commander.getSignIndex();
    size_t caseIndex = commander.getCaseIndex();
    size_t controlIndex = commander.getControlIndex();
    double signNull = commander.getSignNull();
    double infoThreshold = commander.getINFO();
    int nCase = commander.getNCase();
    int nControl = commander.getNControl();
    int nSample = commander.getNSample();
    size_t maxIndex =commander.getMaxIndex();
    //We have got all the required information, now read the file and start processing

    size_t nSkipped = 0, nDuplicate=0, nInvalid=0, nFilter=0, nOverSig = 0;
    std::map<std::string, bool> duplication;

    while(getline(pValue, line)){
        line = usefulTools::trim(line);
        if(!line.empty()){ //Only process lines with information
            std::vector<std::string> token;
            usefulTools::tokenizer(line, " \t", &token);
            if(token.size() > maxIndex){
                std::string rsId = token[rsIndex-1];
                if(duplication.find(rsId) !=duplication.end()){
                    // Check if the same rs ID has been used
                    nDuplicate++;
                }
                else{ //If not duplicated, we start reading in the information
                    std::string chr = token[chrIndex-1];
                    size_t bp = atoi(token[locIndex-1].c_str());
                    std::string refAllele = (refIndex !=0)?token[refIndex-1]:"";
                    std::string altAllele = (altIndex !=0)?token[altIndex-1]:"";

                    bool invalid =false;
                    double imputeScore = 0.0;
                    if(imputeIndex!=0){
                        if(usefulTools::isNumeric(token[imputeIndex-1])) imputeScore = atof(token[imputeIndex-1].c_str());
                        else invalid = true;
                    }
                    int signOfStat = 0;
                    if(signIndex!=0){
                        if(usefulTools::isNumeric(token[signIndex-1])) signOfStat = ((atof(token[signIndex-1].c_str())< signNull)? -1: 1);
                        else invalid = true;
                    }
                    int sizeOfSample = 0;
                    if(qt && sampleIndex !=0){
                        if(usefulTools::isNumeric(token[sampleIndex-1])) sizeOfSample =atoi(token[sampleIndex-1].c_str());
                        else invalid = true;
                    }
                    else if(qt) sizeOfSample = nSample;

                    size_t sizeOfCase=0;
                    if(bt && caseIndex !=0){
                        if(usefulTools::isNumeric(token[caseIndex-1])) sizeOfCase =atoi(token[caseIndex-1].c_str());
                        else invalid = true;
                    }
                    else if(bt) sizeOfCase = nCase;

                    size_t sizeOfControl = 0;
                    if(bt && controlIndex != 0){
                        if(usefulTools::isNumeric(token[controlIndex-1])) sizeOfControl =atoi(token[controlIndex-1].c_str());
                        else invalid = true;
                    }
                    else if(bt) sizeOfControl = nControl;

                    if(!usefulTools::isNumeric(token[sumStatIndex-1])) invalid = true;

                    if(invalid) nInvalid++; //This line is invalid because some of the numerical information are not numeric
                    else if(imputeScore < infoThreshold) nFilter ++;
                    else{
                        double statistic = atof(token[sumStatIndex-1].c_str());
                        // We cannot convert a p-value of 0 to a valid summary statistic
                        // So we need to remove it
                        if(statistic == 0.0 && isP ) nOverSig++; //This is not a safe comparison as double == is always a problem
                        else{
                            snpList.push_back(new Snp(chr,rsId, bp, sizeOfSample, sizeOfCase, sizeOfControl, refAllele, altAllele, statistic, imputeScore, signOfStat));
                            if(isP) snpList.back().computeSummaryStatistc(); //Otherwise, we have already got the required summary statistics
                            duplication[rsId] = true;
                        }
                    }
                }
            }
            else{
                //This mean the current line is malformed
                nSkipped++;
            }
        }

    }
    pValue.close();
    snpList.sort(Snp::sortSnp);

    std::cerr << std::endl;
    fprintf(stderr, "\nP-Value File Information:\n");
    fprintf(stderr, "==================================\n");
    if(nSkipped != 0)
        fprintf(stderr, "%lu line(s) with insufficient columns\n", nSkipped);
    if(nDuplicate!=0)
        fprintf(stderr, "%lu duplicated SNPs\n", nDuplicate);
    if(nInvalid!=0)
        fprintf(stderr, "%lu invalid SNPs\n", nInvalid);
    if(nFilter!=0)
        fprintf(stderr, "%lu SNPs have info score < %f\n", nFilter, infoThreshold);
    if(nOverSig != 0)
        fprintf(stderr, "%lu SNPs have p-value of 0 and are filtered", nOverSig);
    fprintf(stderr, "%lu SNPs are included from the file\n", snpList.size());
    if(snpList.size() ==0) throw std::runtime_error("ERROR: No SNP is provided");

}


bool Snp::sortSnp (const Snp& i, const Snp& j){
    if(i.getChr().compare(j.getChr()) == 0)
		if(i.getLoc() == j.getLoc())
			return i.getRs().compare(j.getRs()) < 0;
		else
			return i.getLoc() < j.getLoc();
	else return (i.getChr().compare(j.getChr()) < 0);
}

