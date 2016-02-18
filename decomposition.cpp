#include "decomposition.h"


Decomposition::Decomposition(size_t thread):m_thread(thread){}

Decomposition::~Decomposition()
{
    //dtor
}

void Decomposition::decompose(Linkage &linkage, std::list<size_t> &snpLoc, std::list<size_t>::iterator startDecompIter, std::list<size_t>::iterator endDecompIter, std::list<size_t>::iterator startCopyIter, std::list<size_t>::iterator endCopyIter, std::list<size_t>::iterator startVarIter, std::list<size_t>::iterator endVarIter, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> &regionList, bool sign, bool start, bool ending){
    // Might want to add a check here to check if the iterators are out of bound

    size_t sizeOfMatrix = std::distance(startDecompIter, endDecompIter);
    size_t copyStartCoordinate = std::distance(snpLoc.begin(), startCopyIter);
    size_t startCoordinate = std::distance(snpLoc.begin(), startDecompIter);
    arma::vec fStat = arma::vec(sizeOfMatrix, arma::fill::zeros);
    arma::vec zStat = arma::vec(sizeOfMatrix, arma::fill::zeros);
    arma::vec nSample = arma::vec(sizeOfMatrix, arma::fill::zeros);
    arma::vec heritResult = arma::vec(sizeOfMatrix, arma::fill::zeros);
    size_t i = 0;
    for(std::list<size_t>::iterator snpLocIter = startDecompIter; snpLocIter != endDecompIter; ++snpLocIter){
        double stat = snpList.at(*(snpLocIter)).getStat();
        int sampleSize = snpList.at(*(snpLocIter)).getSampleSize();
        zStat(i) = stat;
        fStat(i) = (stat*stat-1.0)/((double)sampleSize-2.0+stat*stat); // This is the f-statistic
        nSample(i) = (double)1/(double)sampleSize;
        ++i;
    }

    if(sign){
        arma::mat varResult = arma::mat(sizeOfMatrix,sizeOfMatrix,arma::fill::eye);
        linkage.decompose(startCoordinate,zStat, fStat, nSample, heritResult, varResult);
        /** Get the heritability **/
        std::list<size_t>::iterator iter = startCopyIter;
        for(size_t j = copyStartCoordinate; j < sizeOfMatrix; ++j){
            snpList.at(*(iter)).setHeritability(heritResult(j));
            std::advance(iter,1);
        }
//        std::cerr << "Done copying" << std::endl;
        /** Get the variance **/
        /** This is where the true problem occurs **/
        /** We actually need to know if this is the end **/
        /** Moreover, we need to know about the mid band's location **/
        // Use three for loop, to minimize the conditioning
        // Use iterator for the loop for easy checking
        size_t xCor = 0, yCor = 0;

        xCor = (start)? 0:std::distance(startDecompIter,endVarIter);
//        std::cerr << "Copy first part of variance " << std::endl;


        for(std::list<size_t>::iterator iter=startDecompIter; iter != startVarIter; ++iter){
            // This is the top part
            // For this part, we only take the ending
            size_t currentX = xCor;
            for(std::list<size_t>::iterator innerIter=(start)?startDecompIter:endVarIter; innerIter != endDecompIter; ++innerIter){
                // Now we need to check if both SNPs belongs to the region
                // Only add it to the region if both belongs to the same region
                for(size_t i = 0; i < regionList.size(); ++i){
                    if(snpList.at(*(iter)).flag(i)&&snpList.at(*(innerIter)).flag(i)){
                        regionList[i].addVariance(varResult(currentX,yCor));
                    }
                }
                currentX++;
            }
            yCor++;
        }

        for(std::list<size_t>::iterator iter=startVarIter; iter != endVarIter; ++iter){
            // This is for the mid band
            // For this part, we take everything
            size_t currentX = 0;
            for(std::list<size_t>::iterator innerIter=startDecompIter; innerIter != endDecompIter; ++innerIter){
                // Now we need to check if both SNPs belongs to the region
                // Only add it to the region if both belongs to the same region
                for(size_t i = 0; i < regionList.size(); ++i){
                    if(snpList.at(*(iter)).flag(i)&&snpList.at(*(innerIter)).flag(i)){
                        regionList[i].addVariance(varResult(currentX,yCor));
                    }
                }
                currentX++;
            }
            yCor++;
        }


        std::list<size_t>::iterator copyEnd = (ending)?endDecompIter:startVarIter;
        for(std::list<size_t>::iterator iter=endVarIter; iter != endDecompIter; ++iter){
            // This is the last part
            // For this part, we only take the beginning stuffs
            size_t currentX = 0;
            for(std::list<size_t>::iterator innerIter=startDecompIter; innerIter != copyEnd; ++innerIter){
                for(size_t i = 0; i < regionList.size(); ++i){
                    if(snpList.at(*(iter)).flag(i)&&snpList.at(*(innerIter)).flag(i)){
                        regionList[i].addVariance(varResult(currentX,yCor));
                    }
                }
                currentX++;
            }
            yCor++;
        }
//        std::cerr << "Finished all copying stuff" << std::endl;


    }
    else{
        arma::vec varResult = arma::vec(sizeOfMatrix, arma::fill::zeros);
        linkage.decompose(startCoordinate,fStat,heritResult, varResult);
        // In this case, we can just update the effective number
        /** Get the effective number and heritability **/
        std::list<size_t>::iterator iter = startCopyIter;
        for(size_t j = copyStartCoordinate; j < sizeOfMatrix; ++j){
            snpList.at(*(iter)).setHeritability(heritResult(j));
            snpList.at(*(iter)).setEffective(varResult(j));
            std::advance(iter,1);
        }
    }
}

void Decomposition::run(Linkage &linkage, std::list<size_t> &snpLoc, std::deque<std::list<size_t>::iterator > &boundary, boost::ptr_vector<Snp> &snpList, bool finalizeBuff, bool decomposeAll, bool starting, boost::ptr_vector<Region> &regionList){
    // Here is the meat of the decomposition preprocessing
    // The decomposition is actually done by the linkage class
    // First, get the vectors ready

    // There are a number of situations:
    // 1. Only 1 block is given
    // 2. Only 2 blocks are given
    // 3. At least 3 blocks are given
    bool sign = !(snpList.front().getSign()==0); // Check whether if sign is given

    /** KEY POINT HERE: We don't need to know if it is the last one, because if it isn't the last
     *                  the information will be overwrote in the upcoming round
     */
    if(starting){
        // Check if we decomposeAll, because the situation may differ in that case
        if(decomposeAll){
            // Now this depends on the boundary size
            if(boundary.size() < 4){
                // Then we do the decomposition on everything and update everything
                decompose(linkage, snpLoc,
                          snpLoc.begin(), snpLoc.end(), // Define the decomposition block
                          snpLoc.begin(), snpLoc.end(), // Define the heritability vector copy range
                          snpLoc.begin(), snpLoc.end(), // Define the mid band region of the variance matrix
                          snpList, regionList, sign, starting, true);
            }
            else if(boundary.size() ==4){
                // Then we do the decomposition on the first 3 blocks
                // After that, we decompose the last block

                decompose(linkage, snpLoc,
                          snpLoc.begin(), boundary.back(),  // Define the decomposition block
                          snpLoc.begin(), boundary[2],      // Define the heritability vector copy range
                          boundary[1], boundary[2],      // Define the mid band region of the variance matrix
                          snpList, regionList, sign, starting, false);

                decompose(linkage, snpLoc,
                          boundary[1], snpLoc.end(),    // Define the decomposition block
                          boundary[2], snpLoc.end(),    // Define the heritability vector copy range
                          boundary[2], boundary.back(), // Define the mid band region of the variance matrix
                          snpList, regionList, sign, starting, true);
            }
            else throw std::runtime_error("Undefined condition: check decompose"); //DEBUG message
        }
        else{
            // We will always ignore the last block
            // So decompose everything ignore the last block
            if(boundary.size()-1 < 3 && finalizeBuff){
                //Decompose everything except the last
                decompose(linkage, snpLoc,
                      snpLoc.begin(), boundary.back(),  // Define the decomposition block
                      snpLoc.begin(), boundary.back(),  // Define the heritability vector copy range
                      snpLoc.begin(), boundary.back(),      // Define the mid band region of the variance matrix
                      snpList, regionList, sign, starting, finalizeBuff);
            }
            else if(boundary.size()-1<3) throw std::runtime_error("Must have 4 boundaries if this is the start and not the end");
            else if(boundary.size()==4){
                decompose(linkage, snpLoc,
                      snpLoc.begin(), boundary.back(),  // Define the decomposition block
                      snpLoc.begin(), boundary.back(),  // Define the heritability vector copy range
                      boundary[1], boundary[2],      // Define the mid band region of the variance matrix
                      snpList, regionList, sign, starting, finalizeBuff);
            }
            else throw std::runtime_error("Undefined condition!");
        }
    }
    else{
        // In this case, we know we will at least have 3 blocks
        if(boundary.size() < 3) throw std::runtime_error("For non-starting window(s), there must be at least 3 blocks");
        if(decomposeAll && boundary.size()==3){
            //In this case, decompose everything
//            std::cerr << "First" << std::endl;
            decompose(linkage, snpLoc,
                      snpLoc.begin(), snpLoc.end(), // Define the decomposition block
                      boundary[1], snpLoc.end(),    // Define the heritability vector copy range
                      boundary[1], snpLoc.end(), // Define the mid band region of the variance matrix
                      snpList, regionList, sign, starting, true);
//                      std::cerr << "Fine" << std::endl;
        }
        else if(boundary.size() ==3 && !decomposeAll) throw std::runtime_error("This is impossible for non-start region"); // for my debugging
        else if(boundary.size() ==4){
            // This is most common case, decompose everything except the last block, then continue

//            std::cerr << "second" << std::endl;
            decompose(linkage, snpLoc,
                        snpLoc.begin(), boundary.back(),// Define the decomposition block
                        boundary[1], boundary.back(),   // Define the heritability vector copy range
                        boundary[1], boundary[2],        // Define the mid band region of the variance matrix
                        snpList, regionList, sign, starting, finalizeBuff && !decomposeAll);
//                        std::cerr << "Fine" << std::endl;
            if(decomposeAll){
                // Then we also decompose the last block
//            std::cerr << "end" << std::endl;
                decompose(linkage, snpLoc,
                            boundary[1], snpLoc.end(),  // Define the decomposition block
                            boundary[2], snpLoc.end(),  // Define the heritability vector copy range
                            boundary[2], boundary.back(),
                            snpList, regionList, sign, starting, true);
//                            std::cerr << "Fine" << std::endl;
            }
        }

    }

}
