#include "region.h"


Region::Region(){
    m_names.push_back("With LD");
}

Region::~Region()
{
	//dtor
}

/**
 * Important to note: BED file is 0 base and the end bound is exclusive
 */
void Region::generateRegion(std::string regionList){
    boost::ptr_vector<Interval> padRegion; //For the default region
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
				boost::ptr_vector<Interval> currentRegion;
				m_names.push_back(name);
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
