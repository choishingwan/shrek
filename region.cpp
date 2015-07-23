#include "region.h"


Region::Region(){
    m_names.push_back("With LD");
    m_variance.push_back(0.0);
    m_bufferVariance.push_back(0.0);
}

Region::~Region()
{
	//dtor
}


void Region::Addvariance(double const var, size_t i){
    if(i >= m_variance.size()){
        throw std::out_of_range("Region was out of bound");
    }
    m_variance[i] += var;
}

void Region::CleanBuffer(){
    for(size_t i =0; i < m_bufferVariance.size(); ++i){
        m_bufferVariance[i] = 0;
    }
}


void Region::AddbufferVariance(size_t i, double const var){
    m_bufferVariance.at(i) += var;
}

void Region::SetbufferVariance(double const var, size_t i){
    m_bufferVariance.at(i) = var;

}

void Region::Debuffer(){
    for(size_t i = 0; i < m_bufferVariance.size(); ++i){
        m_variance[i] += m_bufferVariance[i];
        m_bufferVariance[i] = 0;
    }
}

std::string Region::Getname(size_t i) const{
    return m_names.at(i);
}

double Region::Getvariance(double heritability, size_t i, double adjustment) const{
    return adjustment*adjustment*m_variance.at(i);
}

std::string  Region::Getchr(size_t i, size_t j) const{
    return m_intervalList.at(i).at(j)->Getchr();
}

size_t Region::Getstart(size_t i, size_t j) const{
    return m_intervalList.at(i).at(j)->Getstart();
}
size_t Region::Getend(size_t i, size_t j) const{
    return m_intervalList.at(i).at(j)->Getend();
}

size_t Region::GetintervalSize(size_t i) const{
    return m_intervalList.at(i).size();
}
size_t Region::GetnumRegion() const {return m_names.size(); }

void Region::clean(){
    for(size_t i = 0; i < m_intervalList.size(); ++i){
        for(size_t j = 0; j < m_intervalList[i].size(); ++j){
            delete m_intervalList[i][j];
        }
        m_intervalList[i].clear();
    }
    m_intervalList.clear();
}


void Region::generateRegion(std::string regionList){
    std::vector<Interval*> padRegion; //For the default region
    m_intervalList.push_back(padRegion);


    regionList = usefulTools::trim(regionList);
    if(regionList.empty()) return;
    std::vector<std::string> seperateFiles;
    usefulTools::tokenizer(regionList, ", ", &seperateFiles);
    for(size_t i = 0; i < seperateFiles.size(); ++i){
        std::vector<std::string> fileInfo;
        usefulTools::tokenizer(seperateFiles[i], ":", &fileInfo);
        if(fileInfo.size() != 2){
            std::cerr << "Invalid region format, will ignore this region: " << seperateFiles[i] << std::endl;
        }
        else{
            std::string name = fileInfo[0];
            std::string fileName = fileInfo[1];
            std::ifstream regionFile;
            regionFile.open(fileName.c_str());
            if(!regionFile.is_open()){
				std::cerr << "Cannot open region file: " <<fileName << std::endl;
                std::cerr << "Will skip this region" << std::endl;
            }
            else{
				std::vector<Interval* > currentRegion;
				m_names.push_back(name);
                m_variance.push_back(0.0);
                m_bufferVariance.push_back(0.0);
                std::string line;
                while(std::getline(regionFile, line)){
                    line = usefulTools::trim(line);
					std::vector<std::string> bedToken;
                    usefulTools::tokenizer(line, "\t ", &bedToken );
                    if(bedToken.size() > 2){
						currentRegion.push_back(new Interval(bedToken[0], atoi(bedToken[1].c_str()), atoi(bedToken[2].c_str())));
					}
                }
                regionFile.close();
				m_intervalList.push_back(currentRegion);
            }
        }
    }

}
