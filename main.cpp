#include <iostream>
#include <exception>
#include <stdexcept>
#include <ctime>
#include <stdio.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <omp.h>
#include "command.h"
#include "region.h"
#include "snp.h"
#include "genotypefilehandler.h"
#include "snpestimation.h"
//#include "riskprediction.h"

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


extern "C" void openblas_set_num_threads(int); //This is for controlling the multi-thread
extern "C" void goto_set_num_threads(int);

int main(int argc, char *argv[]){
    //Parsing the parameters
    //Start working one by one
    try{
        //We first parse the command line inputs
        Command commander; //Initialize the class (this will set all the default parameters)
        if(!commander.parseCommand(argc, argv)){ //Parameter parsing
            return EXIT_SUCCESS;
        }
        openblas_set_num_threads(commander.getNThread());
        goto_set_num_threads(commander.getNThread());
        //Printing out the input parameters
        //commander.printRunSummary(std::to_string(regionInfo.getNumRegion()));
        time_t now = time(0);
        char* dt = ctime(&now);
        std::cerr << std::endl << "Time: " << dt << std::endl;
        std::cerr << "Command: " << argv[0] << " " << argv[argc-1];
        for(int i = 1; i < argc -1; ++i){
            std::cerr << " " << argv[i];
        }
        std::cerr << std::endl; //This is the simplest way to present the input
        double startTimer = omp_get_wtime();
        //Now we need to parse the region information
        //Parsing the region information
        boost::ptr_vector<Region> regionList;
        Region::generateRegionList(regionList, commander.getRegion());

        // The regionList should now contain all the region information
        // Now start to get the input from the p-value file
        boost::ptr_vector<Snp> snpList;
        // We first get the SNP objects
        Snp::generateSnpList(snpList, commander);
        // Then we build the index of for the SNPs
        std::map<std::string, size_t> snpIndex;
        Snp::generateSnpIndex(snpIndex,snpList,regionList);

        // Next, we need to read the reference panel and generate the required block information
        // Also, we will know what SNPs to be included in the decomposition
        // This is where the first major change occurs
        // 1. We will change to Armadillo, hoping to increase in speed (although the compilation might be a problem)
        // 2. We will also remove ALL the SNPs in perfect LD with each other to improve the conditional number.
        //    Hopefully this will help us to obtain better results.
        GenotypeFileHandler genotypeFileHandler;
        genotypeFileHandler.initialize(commander, snpIndex, snpList);

        if(commander.quantitative() || commander.binary()){
            //Now perform the SNP heritability estimation
            SnpEstimation snpEstimation(commander);
            //The estimate command will basically perform all the analysis except the printing of the results
            snpEstimation.estimate(genotypeFileHandler, snpIndex, snpList, regionList);
            snpEstimation.result(snpList,regionList);
        }
        else{
            fprintf(stderr, "We currently only support SNP heritability estimation in quantitative traits or binary traits\n");
            return EXIT_FAILURE;
        }

        double endTimer = omp_get_wtime();
        fprintf(stderr, "\nCompleted\n");
        fprintf(stderr, "Took %f seconds\n\n", endTimer-startTimer);

    }
    catch( char const* error){
        std::cerr << error << std::endl;
        return EXIT_FAILURE;
    }
    catch(const std::runtime_error& error){
        std::cerr << error.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (std::bad_alloc& ba){
        std::cerr << ba.what() <<std::endl;
        return EXIT_FAILURE;
    }

    return 0;
}
