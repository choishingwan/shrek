#include "riskprediction.h"


RiskPrediction::RiskPrediction(const Command *commander,std::vector<Snp*> *snpList):m_snpList(snpList){
    m_thread = commander->Getthread();
	m_minBlock = commander->GetminBlock();
    m_maxBlock  = commander->GetmaxBlock();
    m_maf = commander->Getmaf();
    m_ldCorrection = commander->ldCorrect();
    m_keep = commander->keep();
    m_maxBlockSet = commander->maxBlockSet();
    m_genotypeFilePrefix = commander->GetgenotypeFile();
}


RiskPrediction::~RiskPrediction()
{
    //dtor
}


RiskPrediction::checkGenotype(){
    std::string bimFileName = m_genotypeFilePrefix;
    bimFileName.append(".bim");
    std::ifstream bimFile;
    bimFile.open(bimFileName.c_str());
    if(!bimFile.is_open()){
        std::string message = "Cannot open bim file: ";
        message.append(bimFileName);
        throw message;
    }
    std::string line;
    while(std::getline(bimFile, line)){
        line = usefulTools::trim(line);
        std::vector<std::string> token;
        usefulTools::tokenizer(line, "\t ", &token);
        if(!line.empty() && token.size() >= 6 ){
            std::string chr = token[0];
            std::string rsId = token[1];
            size_t loc = atoi(token[3].c_str());
            std::string ref = token[4];
            std::string alt = token[5];

        }
    }
    bimFile.close();
}
