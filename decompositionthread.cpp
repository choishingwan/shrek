#include "decompositionthread.h"

std::mutex DecompositionThread::decomposeMtx;
Eigen::MatrixXd DecompositionThread::checking;


DecompositionThread::DecompositionThread(size_t start, size_t length, Eigen::VectorXd const * const betaEstimate, Eigen::VectorXd const * const sqrtChiSq, Linkage *linkage, std::deque<size_t>  *snpLoc, std::vector<Snp*> *snpList, bool chrStart, bool lastOfBlock, bool secondLastOfBlock,  Region *regionInfo):m_start(start), m_length(length), m_betaEstimate(betaEstimate), m_sqrtChiSq(sqrtChiSq),m_linkage(linkage), m_snpLoc(snpLoc), m_snpList(snpList), m_chrStart(chrStart), m_lastOfBlock(lastOfBlock), m_secondLastOfBlock(secondLastOfBlock), m_regionInfo(regionInfo){}

DecompositionThread::~DecompositionThread()
{}


void *DecompositionThread::ThreadProcesser(void *in){
    ((DecompositionThread *) in)->solve();
	return nullptr;
}

void DecompositionThread::solve(){
    //Here is where I need to change stuff
    size_t betaLength = (*m_betaEstimate).rows();
    size_t processLength = m_length;
    if(m_lastOfBlock) processLength=betaLength-m_start;
    Eigen::MatrixXd variance = Eigen::MatrixXd::Zero(processLength, processLength);
    Eigen::MatrixXd additionVariance = Eigen::MatrixXd::Zero(processLength, processLength);
	Eigen::VectorXd result;
    result = m_linkage->solve(m_start, processLength, m_betaEstimate, m_sqrtChiSq, &variance, &additionVariance, Snp::GetmaxSampleSize());;
    double multiplier = 2.0;
	size_t copyStart = m_length/3;
	size_t copyEnd = m_length/3;
	if(m_chrStart && m_start == 0){
		copyStart = 0;
        copyEnd += m_length/3;
        multiplier=1.0;
	}
	if(m_lastOfBlock) multiplier = 1.0; //So that we can account for the varying end block size.


    std::vector<double> regionVariance(m_regionInfo->GetnumRegion(), 0.0);
    std::vector<double> regionAdditionVariance(m_regionInfo->GetnumRegion(), 0.0);
    std::vector<double> regionBufferVariance(m_regionInfo->GetnumRegion(), 0.0);
    std::vector<double> regionBufferAdditionVariance(m_regionInfo->GetnumRegion(), 0.0);


    for(size_t i = 0; i < m_length/3; ++i){
        for(size_t j = copyEnd+copyStart; j < processLength; ++j){
            double covariance = variance(i,j);
			double additionCovariance = additionVariance(i,j);
			for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                if((*m_snpList)[(*m_snpLoc)[m_start+i]]->GetFlag(regionIndex)&&(*m_snpList)[(*m_snpLoc)[m_start+j]]->GetFlag(regionIndex)){
					if(m_lastOfBlock || m_secondLastOfBlock){
                        //For the last 2 blocks, it is therefore perfect LD and the whole thing should
                        //be stored within the buffer such that only if we have reached the end
                        //that we should put these buffer back into the region variance
                        regionBufferVariance[regionIndex]+= multiplier*covariance;
                        regionBufferAdditionVariance[regionIndex]+= multiplier*additionCovariance;
					}
					else{
                        regionVariance[regionIndex]+= multiplier*covariance;
                        regionAdditionVariance[regionIndex]+= multiplier*additionCovariance;
					}
                }
			}

        }
    }

    for(size_t i = copyStart; i < copyStart+copyEnd; ++i){
		(*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability(result(i));
		(*m_snpList)[(*m_snpLoc)[m_start+i]]->Setvariance(variance(i,i)); //The diagonal of the matrix contain the per snp variance
		(*m_snpList)[(*m_snpLoc)[m_start+i]]->SetadditionVariance(additionVariance(i,i)); //The diagonal of the matrix contain the per snp variance
        for(size_t j = 0; j < processLength; ++j){ //For this snp, go through all the snp partners

			double covariance = variance(i,j);
			double additionCovariance = additionVariance(i,j);
			//checking(m_start+i,m_start+j) += covariance; //DEBUG
			for(size_t regionIndex = 0; regionIndex < regionVariance.size(); ++regionIndex){
                if((*m_snpList)[(*m_snpLoc)[m_start+i]]->GetFlag(regionIndex)&&
                   (*m_snpList)[(*m_snpLoc)[m_start+j]]->GetFlag(regionIndex)){
                    if(m_lastOfBlock || m_secondLastOfBlock){
                        //same as above
                        regionBufferVariance[regionIndex]+= covariance;
                        regionBufferAdditionVariance[regionIndex]+= additionCovariance;
					}
					else{
                        regionVariance[regionIndex]+= covariance;
                        regionAdditionVariance[regionIndex]+= additionCovariance;
					}
                }
			}
        }
	}
// TODO (swchoi#1#): Something is still wrong with the registration of the variance. Need to find out the reason and tackle it accordingly

    /** Now for the rest, this should only be use when it is the last block. If
     *  this is the last block, then instead of putting the bottom portion into
     *  the actual variance and addition variance vector, we put them to the buffer
     *  which will only be sent to the actual vector when we encountered a new
     *  chromosome or the end of the process. Therefore avoiding the double counting
     *  of those portion of variance.
     */
    if(m_lastOfBlock){
        //This is only applicable for the last block but not the second last block
        for(size_t i = copyStart+copyEnd; i < processLength; ++i){
            /** Doesn't matter with these three as if we are indeed not the last
             *  block of the chromosome, they will be re-wrote
             */
            (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setheritability(result(i));
            (*m_snpList)[(*m_snpLoc)[m_start+i]]->Setvariance(variance(i,i));
            (*m_snpList)[(*m_snpLoc)[m_start+i]]->SetadditionVariance(additionVariance(i,i));
            for(size_t j =  0; j < processLength; ++j){
                double covariance = variance(i,j);
                double additionCovariance = additionVariance(i,j);
                //checking(m_start+i,m_start+j) += covariance; //DEBUG
                for(size_t regionIndex = 0; regionIndex < m_regionInfo->GetnumRegion(); ++regionIndex){
                    if((*m_snpList)[(*m_snpLoc)[m_start+i]]->GetFlag(regionIndex)&&
                       (*m_snpList)[(*m_snpLoc)[m_start+j]]->GetFlag(regionIndex)){
                        regionBufferVariance[regionIndex]+= covariance;
                        regionBufferAdditionVariance[regionIndex]+= additionCovariance;
                    }
                }
            }
        }
    }

	decomposeMtx.lock();
		for(size_t i = 0; i < m_regionInfo->GetnumRegion(); ++i){
            m_regionInfo->Addvariance(regionVariance[i], i);
            m_regionInfo->AddadditionVariance(regionAdditionVariance[i], i);
            m_regionInfo->SetbufferVariance(regionBufferVariance[i],i);
            m_regionInfo->SetbufferAdditionVariance(regionBufferAdditionVariance[i], i);
		}
	decomposeMtx.unlock();

}
