#include <iostream>
#include <exception>
#include <stdexcept>
#include <boost/ptr_container/ptr_vector.hpp>
#include "command.h"
#include "region.h"
#include "snp.h"
#include "genotypefilehandler.h"
#include "snpestimation.h"
#include "riskprediction.h"

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
*   @version 0.02
*/




int main(int argc, char *argv[]){
    //Parsing the parameters
    //Start working one by one
    try{
        Command *commander = nullptr;
        commander=new Command();
        commander->initialize(argc, argv);
        Region *regionInfo = nullptr;
        regionInfo = new Region();
        regionInfo->generateRegion(commander->getRegion());
        commander->printRunSummary(std::to_string(regionInfo->getNumRegion()));
        boost::ptr_vector<Snp> snpList;
        std::map<std::string, size_t> snpIndex;
        Snp::generateSnpList(snpList, commander);
        Snp::generateSnpIndex(snpIndex, snpList, commander, regionInfo);
        boost::ptr_vector<Interval> blockInfo; //We use this to store all block information, should be useful for both risk and not risk stuff
        if(commander->quantitative() || commander->caseControl()){
            GenotypeFileHandler *genotypeFileHandler = nullptr;
            genotypeFileHandler = new GenotypeFileHandler();
            genotypeFileHandler->initialize(commander, snpIndex, snpList, blockInfo);
            //Now everything is prepared, we can start the SNP heritability estimation
            SnpEstimation *snpEstimation = nullptr;
            snpEstimation = new SnpEstimation();
            snpEstimation->Estimate(genotypeFileHandler, snpIndex, snpList, regionInfo, commander, blockInfo);

        }
        else if(commander->diRisk() || commander->conRisk()){

        }
    }
    catch(const std::runtime_error& error){
        std::cerr << error.what() << std::endl;
    }
    catch (std::bad_alloc& ba){
        std::cerr << ba.what() <<std::endl;
    }


/*

    if(commander->risk()){
        RiskPrediction *riskPrediction = nullptr;
        riskPrediction = new RiskPrediction(commander,&snpList);
        try{
            riskPrediction->checkGenotype();
            riskPrediction->run();
            riskPrediction->result();
            return EXIT_SUCCESS;
        }
        catch(const char *e){
            std::cerr << e << std::endl;
            delete commander;
            regionInfo->clean();
            delete regionInfo;
            Snp::cleanSnp(snpList);
            delete riskPrediction;
            return EXIT_FAILURE;
        }
    }
    else{
        try{
        */
            /** Generate the Snp Index */
            /*
            if(commander->quantitative()){
                Snp::generateSnpIndex(snpIndex, snpList, regionInfo, commander->isPvalue(), commander->GetextremeAdjust());
            }
            else if(commander->caseControl()){
                Snp::generateSnpIndex(snpIndex, snpList,commander->GetcaseSize(), commander->GetcontrolSize(), commander->Getprevalence(), regionInfo, commander->isPvalue());
                */
                /** For case control study, set the liability adjustment */
                /*
                Snp::Setadjustment(commander->Getprevalence(), commander->GetcaseSize(), commander->GetcontrolSize());
            }
            regionInfo->clean();
        }
        catch(const char *e){
            std::cerr << e << std::endl;
            delete commander;
            regionInfo->clean();
            delete regionInfo;
            Snp::cleanSnp(snpList);
            return EXIT_FAILURE;
        }
        if(!commander->GetdirectionFile().empty()){
            try{
            */
                /** If there is direction information, use it */
                /*
                Snp::addDirection(snpIndex, snpList, commander->GetdirectionFile());
            }
            catch(const char *e){
                std::cerr << e << std::endl;
                delete commander;
                delete regionInfo;
                Snp::cleanSnp(snpList);
                return EXIT_FAILURE;
            }
        }
        //From now on, we are only allow to iterate through snpList through snpIndex
        GenotypeFileHandler *genotypeFileHandler = nullptr;
        genotypeFileHandler = new GenotypeFileHandler(commander->GetldFilePrefix(), commander->Getthread(), commander->GetoutputPrefix());
        try{
            genotypeFileHandler->initialize(snpIndex, &snpList, commander->validate(), commander->maxBlockSet(), commander->GetmaxBlock(), commander->GetminBlock(),commander->Getmaf());
        }
        catch (const char *e) {
            std::cerr << "Exception encountered when opening genotype files" << std::endl;
            std::cerr << e << std::endl;
            delete commander;
            delete regionInfo;
            Snp::cleanSnp(snpList);
            delete genotypeFileHandler;
            return EXIT_FAILURE;
        }

        SnpEstimation *snpEstimation = nullptr;
        snpEstimation = new SnpEstimation(genotypeFileHandler, &snpIndex, &snpList, commander->Getthread(), commander->Getmaf(), commander->ldCorrect(), regionInfo);

        try{
            snpEstimation->Estimate();
        }
        catch(const char *e){
            std::cerr << e << std::endl;
            delete commander;
            delete snpEstimation;
            Snp::cleanSnp(snpList);
            delete regionInfo;
            return EXIT_FAILURE;
        }

        snpEstimation->Getresult(commander->GetoutputPrefix());
        delete snpEstimation;
        delete genotypeFileHandler;
    }
    //Cleaning section
	Snp::cleanSnp(snpList);
    delete commander;
    regionInfo->clean();
    delete regionInfo;
    */
    return 0;
}
