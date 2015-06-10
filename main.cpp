#include <iostream>
#include "command.h"
#include "region.h"
#include "snp.h"
#include "genotypefilehandler.h"
#include "snpestimation.h"

/** @mainpage SHREK: Snp HeRitability Estimate Kit
*   @par Description:
*	This is the heritability estimation kit. This programme uses the test-statistic from GWAS studies
*	to estimate the heritability of a certain trait.
*	Details of the method, limitations and other information can be found in the corresponding read me.
*	@author Sam S.W. Choi
*	@author Johnny S.H. Kwan
*	@author Henry C.M. Leung
*	@author Pak C. Sham
*	@author Tim S.H. Mak
*	@author D. Campbell
*	@author M.X. Li
*   @version 0.01
*/

// TODO (swchoi#1#): Need to find a way to output the perfect SNP pairs
// TODO (swchoi#1#): Update help message in the Command Class ...




int main(int argc, char *argv[]){
    //Parsing the parameters
    Command *commander = new Command();
    try{
        /** Parse the parameters using the command handler */
		commander->initialize(argc, argv);
    }
    catch(const char* e){
		std::cerr << e << std::endl;
		delete commander;
		return EXIT_FAILURE;
    }
    catch(const int e){
		delete commander;
        return EXIT_SUCCESS;
    }
    Region *regionInfo = new Region();
    try{
        regionInfo->generateRegion(commander->GetregionList());
    }
    catch(const char *e){
        std::cerr << e << std::endl;
        regionInfo->clean();
        delete commander;
        delete regionInfo;
        return EXIT_FAILURE;
    }
    commander->printRunSummary(std::to_string(regionInfo->GetnumRegion()));
    std::vector<Snp*> snpList;
    SnpIndex *snpIndex = new SnpIndex();
    try{
        Snp::generateSnpList(snpList, commander);
    }
    catch (const std::ifstream::failure e) {
        std::cerr << "Exception encountered when opening file" << std::endl;
        std::cerr << "Please check all your input file are ok" << std::endl;
        delete commander;
        regionInfo->clean();
        delete regionInfo;
        delete snpIndex;
        Snp::cleanSnp(snpList);
        return EXIT_FAILURE;
    }
    catch(const char *e){
        std::cerr << e << std::endl;
        delete commander;
        regionInfo->clean();
        delete regionInfo;
        delete snpIndex;
        Snp::cleanSnp(snpList);
        return EXIT_FAILURE;
    }
    try{
        if(commander->quantitative()){
            Snp::generateSnpIndex(snpIndex, snpList, regionInfo, commander->isPvalue(), commander->GetextremeAdjust());
        }
        else if(commander->caseControl()){
            Snp::generateSnpIndex(snpIndex, snpList,commander->GetcaseSize(), commander->GetcontrolSize(), commander->Getprevalence(), regionInfo, commander->isPvalue());
        }
        regionInfo->clean();
    }
    catch(const char *e){
        std::cerr << e << std::endl;
        delete commander;
        regionInfo->clean();
        delete regionInfo;
        delete snpIndex;
        Snp::cleanSnp(snpList);
        return EXIT_FAILURE;
    }

	//From now on, we are only allow to iterate through snpList through snpIndex
    GenotypeFileHandler *genotypeFileHandler = new GenotypeFileHandler(commander->GetldFilePrefix(), commander->Getthread());
    try{
        genotypeFileHandler->initialize(snpIndex, &snpList, commander->validate(), commander->maxBlockSet(), commander->GetmaxBlock(), commander->GetminBlock());
    }
    catch (const std::ifstream::failure e) {
        std::cerr << "Exception encountered when opening genotype files" << std::endl;
        delete commander;
        delete regionInfo;
        delete snpIndex;
        Snp::cleanSnp(snpList);
        delete genotypeFileHandler;
        return EXIT_FAILURE;
    }

	SnpEstimation *snpEstimation = new SnpEstimation(genotypeFileHandler, snpIndex, &snpList, commander->Getthread(), commander->Getmaf(), commander->ldCorrect(), regionInfo);
    try{
        snpEstimation->Estimate();
    }
    catch(const char *e){
        std::cerr << e << std::endl;
        delete commander;
        delete snpEstimation;
        delete snpIndex;
        Snp::cleanSnp(snpList);
        delete regionInfo;
        return EXIT_FAILURE;
    }

    snpEstimation->Getresult(commander->GetoutputPrefix());

    //Cleaning section
	Snp::cleanSnp(snpList);
	delete snpEstimation;
	delete snpIndex;
	delete genotypeFileHandler;
    delete commander;
    regionInfo->clean();
    delete regionInfo;

    return 0;
}
