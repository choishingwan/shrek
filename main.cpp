#include <iostream>
#include <exception>
#include <cstdlib>
#include <stdexcept>
#include <ctime>
#include <stdio.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <omp.h>
#include <execinfo.h>
#include <signal.h>
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

void handler(int sig) {
  void *array[10];
  size_t size;
  // get void*'s for all entries on the stack
  size = backtrace(array, 10);
  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]){
    //Parsing the parameters
    //Start working one by one
    signal(SIGSEGV, handler);
    try{
        //We first parse the command line inputs
        Command commander; //Initialize the class (this will set all the default parameters)
        if(!commander.parseCommand(argc, argv)) return EXIT_SUCCESS; //If return false, then it means we don't need to do anything else
        // Here we hope to control for the threads used by openBLAS
        openblas_set_num_threads(commander.getNThread());
        goto_set_num_threads(commander.getNThread());
        time_t now = time(0);
        char* dt = ctime(&now);
        fprintf(stderr, "\nTime: %s\n", dt);
        fprintf(stderr, "Command: %s %s", argv[0], argv[argc-1]);
        for(int i = 1; i < argc -1; ++i) fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n");
        double startTimer = omp_get_wtime(); //Get the start time of the programme
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
        // In this new algorithm, the snpIndex is used as a quick search for the genotypeFileHandler
        std::map<std::string, size_t> snpIndex;
        Snp::generateSnpIndex(snpIndex,snpList,regionList);

        // Next, we need to read the reference panel
        GenotypeFileHandler genotypeFileHandler;
        // Now we will initialize the handler (e.g. check the bim and fam file)
        genotypeFileHandler.initialize(commander, snpIndex, snpList);
        //We no longer need the snpIndex;
        snpIndex.clear();

        if(commander.quantitative() || commander.binary()){
            //Now perform the SNP heritability estimation
            SnpEstimation snpEstimation(commander);
            //The estimate command will basically perform all the analysis except the printing of the results
            snpEstimation.estimate(genotypeFileHandler, snpList, regionList);
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
    catch( char const* error){ // Need to learn better error handling of this try catch
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
