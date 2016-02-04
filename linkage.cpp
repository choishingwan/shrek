#include "linkage.h"

std::mutex Linkage::linkageMtx;

Linkage::Linkage(size_t thread):m_thread(thread){}
Linkage::~Linkage(){}

void Linkage::computeLd(const boost::ptr_list<Genotype> &genotype, const std::list<size_t> &snpLoc, size_t index, size_t horiStart, size_t endOfBlock, boost::ptr_vector<Snp> &snpList, const bool correction, std::vector<size_t> &perfectLd, std::vector<size_t> &calculatedR, std::vector<size_t> &calculatedR2){
    // From horiStart to endOfBlock, calculate LD and record down all the perfect LD's index
    std::vector<size_t> buffPerfect; // This is for a buffer, such that we only need to mtx once in the whole function
    // mtx can take up a lot of time if we repeat it over and over again
    assert(endOfBlock<=genotype.size() && "Should not calculate LD for more SNPs than there is within the block");
    boost::ptr_list<Genotype>::const_iterator  iter = genotype.begin(), endIter = genotype.begin(), currentIter = genotype.begin();
    std::advance(iter, horiStart);
    std::advance(endIter, endOfBlock+1); // We use this as a bound
    std::advance(currentIter, index);
    //linkageMtx.lock();
        //fprintf(stderr, "Check: %lu, %lu, %lu, %lu\n", index, horiStart, endOfBlock, genotype.size());
    //linkageMtx.unlock();
    size_t currentIndex = horiStart;
    for(; iter != endIter; ++iter){
        double r = 0.0;
        double r2 = 0.0;
        (*currentIter).GetbothR((*iter),correction,r, r2);
        // Now here is where we need to think think


        // Also update the index
        currentIndex++;


    }
}


// This function is only responsible for distributing the jobs to individual threads
// it should however, also be responsible for updating the genotype, snpLoc and boundary
// after the removal of perfect LD
void Linkage::startComputeLd(boost::ptr_list<Genotype> &genotype, std::list<size_t> &snpLoc, size_t index, size_t horiStart, size_t endOfBlock, boost::ptr_vector<Snp> &snpList, const bool correction){
    // This will construct the LD of this particular SNP
    // Actually, the following is only applicable for the last block
    size_t startOfBlock = (index > horiStart)? index: horiStart;

    size_t numSnp = endOfBlock-startOfBlock;
    if(numSnp == 0) return; // nothing to calculate
    std::vector<double> calculatedR2(numSnp, 0.0); //Initialize the result vector
    std::vector<double> calculatedR(numSnp, 0.0); //Initialize the result vector
    std::vector<size_t> perfectLd; // the vector containing the location of SNPs to be removed
    std::vector<std::thread> threadStore;
    size_t currentStart = horiStart;
    if(numSnp < m_thread){
        // We can use each thread to calculate the number
        for(size_t i = 0; i < numSnp; ++i){
            // Put each SNP into a thread for the calculation
            threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), startOfBlock, currentStart, currentStart+1, std::ref(snpList), std::cref(correction), std::ref(perfectLd), std::ref(calculatedR), std::ref(calculatedR2)));
            currentStart++;
        }
    }
    else{
        //We have more SNPs than threads, so better divide them into groups
        int remainder = numSnp % m_thread;
        size_t block = numSnp / m_thread;
        for(size_t i = 0; i < m_thread; ++i){
            if(remainder> 0){
                threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), startOfBlock, currentStart, currentStart+block+1, std::ref(snpList), std::cref(correction), std::ref(perfectLd), std::ref(calculatedR), std::ref(calculatedR2)));
                currentStart+= block+1;
                --remainder;
            }
            else{
                threadStore.push_back(std::thread(&Linkage::computeLd, this, std::cref(genotype), std::cref(snpLoc), startOfBlock, currentStart, currentStart+block, std::ref(snpList), std::cref(correction), std::ref(perfectLd), std::ref(calculatedR), std::ref(calculatedR2)));
                currentStart+= block;
            }
        }
    }
    for (size_t i = 0; i < threadStore.size(); ++i) {
        threadStore[i].join();
    }
    threadStore.clear();

    // Now remove the perfect LD stuff here
    std::sort(perfectLd.begin(), perfectLd.end());

}

void Linkage::construct(boost::ptr_list<Genotype> &genotype, std::list<size_t> &snpLoc, std::deque<std::list<size_t>::iterator > &boundary, boost::ptr_vector<Snp> &snpList, const bool correction){
    // This is yet another complicated algorithm
    // The aim of this function is to construct the required LD matrix
    // (Think about it, maybe I will shift the perfect LD removal to the other function)
    if(genotype.empty()){
        throw std::runtime_error("Cannot build LD without genotypes");
	}
    // If this is the first of its kind, then we need to initialize the matrix
    bool start = false;
	if(m_linkage.n_cols==0){
        start = true;
        m_linkage = arma::mat(genotype.size(), genotype.size(),arma::fill::eye);
        m_linkageSqrt = arma::mat(genotype.size(), genotype.size(), arma::fill::eye);
    }
    else{ // Otherwise, just make sure the matrix is of the right size
        m_linkage.resize(genotype.size(), genotype.size());
        m_linkageSqrt.resize(genotype.size(), genotype.size());
    }
    // We don't need the buff information
    // because for this function, we should build all the matrix anyway
    // Now if we are not working on the new stuff, we assume the LD matrix has been clean
    // of the unwanted materials
    // So we only need to construct the matrix here
    // When we construct the LD matrix, we do it per SNP level but still use multi-threading
    // then we skip any SNP that are found to be in perfect LD with others
    // This will allow us to avoid the need of removing rols (only column)
    // which should be fast considering armadillo is column major
    fprintf(stderr, "Start of SNP loc: %lu\n", snpLoc.front());
    fprintf(stderr, "size of boundary: %lu\n", boundary.size());
    fprintf(stderr, "Boundary back: %lu\n", (*(boundary.back())));
    fprintf(stderr, "SNP loc end %lu\n", snpLoc.back());

    if(boundary.size() < 3){
        //Very small input, so we will do the whole part together
        for(size_t i = 0; i < genotype.size(); ++i){
            //Just construct the whole LD matrix
            startComputeLd(genotype, snpLoc, i, i,genotype.size(), snpList, correction);
        }
    }
    else{
        if(start){
            for(size_t i = 0; i < (*(boundary[1])); ++i){
                // Only the first block has some different start and end
                startComputeLd(genotype, snpLoc, i,i, (*(boundary.back())), snpList, correction);
            }
        }
        // All the rest are the same no matter if this is the start or not
        for(size_t i = (*(boundary[1])); i < genotype.size(); ++i){
            if(start)
                startComputeLd(genotype, snpLoc, i, i,genotype.size(), snpList, correction);
            else
                startComputeLd(genotype, snpLoc, i, (*(boundary.back())),genotype.size(), snpList, correction);
        }
    }

}



