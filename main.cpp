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
        //We first parse the command line inputs
        Command commander; //Initialize the class (this will set all the default parameters)
        commander.initialize(argc, argv); //Parameter parsing

        Region regionInfo;
        regionInfo.generateRegion(commander.getRegion());
        commander.printRunSummary(std::to_string(regionInfo.getNumRegion()));
        boost::ptr_vector<Snp> snpList;
        std::map<std::string, size_t> snpIndex;
        Snp::generateSnpList(snpList, commander);
        std::vector<int> genoInclusion; //Bad design here, but don't bother re-writing. useless for heritability estimation
        /**
         * This function will update the snpIndex accordingly and will not include SNPs that are not presented
         * in the genotype file.
         */
        Snp::generateSnpIndex(snpIndex, snpList, commander, regionInfo, genoInclusion);
        boost::ptr_vector<Interval> blockInfo; //We use this to store all block information, should be useful for both risk and not risk stuff
        GenotypeFileHandler genotypeFileHandler;
        genotypeFileHandler.initialize(commander, snpIndex, snpList, blockInfo);
        if(commander.quantitative() || commander.caseControl()){
            //Now everything is prepared, we can start the SNP heritability estimation
            SnpEstimation snpEstimation;
            snpEstimation.Estimate(genotypeFileHandler, snpIndex, snpList, regionInfo, commander, blockInfo);
            snpEstimation.getResult(commander, regionInfo, snpIndex,snpList);
        }
        else if(commander.diRisk() || commander.conRisk()){
            SnpEstimation snpPrediction;
            snpPrediction.Predict(genotypeFileHandler, snpIndex, snpList, regionInfo, commander, blockInfo,genoInclusion);

        }
    }
    catch( char const* error){
        std::cerr << error << std::endl;
    }
    catch(const std::runtime_error& error){
        std::cerr << error.what() << std::endl;
    }
    catch (std::bad_alloc& ba){
        std::cerr << ba.what() <<std::endl;
    }


    return 0;
}
