#include "region.h"

Region::Region(std::string chr, size_t start, size_t end):m_chr(chr), m_start(start), m_end(end){}
std::vector<std::string> Region::regionNames;

std::string Region::Getchr() const { return m_chr; }
size_t Region::Getstart() const { return m_start; }
size_t Region::Getend() const { return m_end; }

Region::~Region()
{
	//dtor
}

void Region::generateRegion(std::vector<std::vector<Region*> > &regionOut, std::string regionList){
	regionList = usefulTools::trim(regionList);
    if(regionList.empty()) return;
    std::vector<std::string> seperateFiles;
    usefulTools::tokenizer(regionList, ",", &seperateFiles);
    if(seperateFiles.size() == 0) return;
    size_t seperateFileSize =seperateFiles.size();
    for(size_t i = 0; i < seperateFileSize; ++i){
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
				std::vector<Region*> currentRegion;
				Region::regionNames.push_back(name);
                std::string line;
                while(std::getline(regionFile, line)){
                    line = usefulTools::trim(line);
					std::vector<std::string> bedToken;
                    usefulTools::tokenizer(line, "\t ", &bedToken );
                    if(bedToken.size() > 2){
						currentRegion.push_back(new Region(bedToken[0], atoi(bedToken[1].c_str()), atoi(bedToken[2].c_str())));
					}
                }
                regionFile.close();
				regionOut.push_back(currentRegion);
            }

        }
    }
}


void Region::cleanRegion(std::vector<std::vector<Region*> > &regionList){
    for(unsigned int i = 0; i < regionList.size(); ++i){
        for(unsigned int j = 0; j < regionList[i].size(); ++j){
            delete regionList[i][j];
        }
        regionList[i].clear();
    }
    regionList.clear();
}
