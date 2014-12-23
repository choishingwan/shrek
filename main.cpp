#include <iostream>
#include "command.h"
#include "region.h"
#include "snp.h"
#include "genotypefilehandler.h"
#include "snpestimation.h"

int main(int argc, char *argv[]){
	Command *commander = new Command(argc, argv);
    std::vector<std::vector<Region*> > regionList;
	Region::generateRegion(regionList, commander->GetregionList());
	commander->printRunSummary(std::to_string(regionList.size()));
	std::vector<Snp*> snpList;
    SnpIndex *snpIndex = new SnpIndex();

	if(commander->quantitative()){
		Snp::generateSnpList(snpList, commander->GetpValueFileName(), commander->GettIndex(), commander->GetsampleSize(), commander->GetrsIndex(), commander->GetbpIndex(), commander->GetchrIndex(), commander->GetsampleSizeIndex(), commander->provideSampleSize());
		Snp::generateSnpIndex(snpIndex, snpList, regionList, commander->isPvalue() );
	}
	else{
		Snp::generateSnpList(snpList, commander->GetpValueFileName(), commander->GetcIndex(), commander->GetsampleSize(), commander->GetrsIndex(), commander->GetbpIndex(), commander->GetchrIndex(), commander->GetsampleSizeIndex(), commander->provideSampleSize());
		 Snp::generateSnpIndex(snpIndex, snpList,commander->GetcaseSize(), commander->GetcontrolSize(), commander->Getprevalence(), regionList, commander->isPvalue() );
	}
	//From now on, we are only allow to iterate through snpList through snpIndex
	GenotypeFileHandler *genotypeFileHandler = new GenotypeFileHandler(commander->GetgenotypeFilePrefix(), snpIndex, snpList);
	SnpEstimation *snpEstimation = new SnpEstimation(genotypeFileHandler, snpIndex, &snpList);

    //Cleaning section
	Snp::cleanSnp(snpList);
	delete snpEstimation;
	delete snpIndex;
	delete genotypeFileHandler;
    delete commander;
    Region::cleanRegion(regionList);
    return 0;
}
