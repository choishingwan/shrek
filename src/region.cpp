#include "region.h"

Region::Region(const std::string name):m_name(name){}

Region::Region(const std::string name,const std::string bedFileName):m_name(name){
    std::ifstream regionFile;
    regionFile.open(bedFileName);
    // This should not happen because we have already checked the regionFile exists before getting into this function
    assert(regionFile.is_open() && "Cannot open region file");
    std::string line;
    while(std::getline(regionFile, line)){
        line= usefulTools::trim(line);
        std::vector<std::string> token;
        usefulTools::tokenizer(line, "\t", &token);
        if(token.size() > 2){
            // The format is chr, start, end
            m_intervalList.push_back(new Interval(token[0], atoi(token[1].c_str()), atoi(token[2].c_str())));
        }
    }
    regionFile.close();
}

Region::~Region(){}


void Region::generateRegionList(boost::ptr_vector<Region> &regionList, const std::string regionInput){
    // At the end, the regionList should contain all the regions
    // The first region contains all the SNPs
    regionList.push_back(new Region("Base"));

    // Check whether if there is any regions
    std::string regionTrim = usefulTools::trim(regionInput);
    if(regionTrim.empty()) return; // if the region input is empty, there is no need to work on it

    // Now get the multiple regions, each should be separated by ','
    std::vector<std::string> seperateFiles;
    usefulTools::tokenizer(regionTrim, ", ", &seperateFiles);
    for(size_t i=0; i < seperateFiles.size(); ++i){
        // For each of the region, extract the region name and form the region
        std::vector<std::string> fileInfo;
        usefulTools::tokenizer(seperateFiles[i], ":", &fileInfo);
        if(fileInfo.size() != 2){
            // The format of the correct region input should always be name:region
            fprintf(stderr,"Invalid region format, will ignore this region: %s\n",seperateFiles[i].c_str());
        }
        else{
            std::string name = fileInfo[0];
            std::string fileName = fileInfo[1];
            if(usefulTools::fileExists(fileName)){
                // Generate the region
                regionList.push_back(new Region(name, fileName));
            }
            else{
                fprintf(stderr, "Cannot open region bed file: %s\n", fileName.c_str());
                fprintf(stderr, "Region skipped\n");
            }

        }
    }
}
