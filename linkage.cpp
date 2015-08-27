#include "linkage.h"
//Eigen::MatrixXd Linkage::m_testing = Eigen::MatrixXd::Zero(3000,3000);

Linkage::Linkage(){
    m_perfectLd =std::vector<size_t>();
    m_snpLoc = nullptr;
    m_snpList = nullptr;
    m_thread = 1;
}

void Linkage::setSnpList(std::vector<Snp*> *snpList){
    m_snpList = snpList;
}
void Linkage::setSnpLoc(std::deque<size_t> *snpLoc){
    m_snpLoc = snpLoc;
}
void Linkage::setThread(size_t thread){
    m_thread = thread;
}

Linkage::~Linkage()
{
	//dtor
}

size_t Linkage::rows() const { return m_linkage.rows(); }
size_t Linkage::cols() const { return m_linkage.cols(); }
Eigen::MatrixXd Linkage::block(size_t blockStart, size_t lengthOfBlock){ return m_linkage.block(blockStart, blockStart, lengthOfBlock, lengthOfBlock); }
Eigen::MatrixXd Linkage::blockSqrt(size_t blockStart, size_t lengthOfBlock){ return m_linkageSqrt.block(blockStart, blockStart, lengthOfBlock, lengthOfBlock); }


void Linkage::triangularThread( const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype){
    //Make the thread region
    std::vector<LinkageThread*> garbageCollection;
    size_t maxThread = (endBlock-startBlock)/ m_thread;
    if(maxThread >=1) maxThread = m_thread;
    else maxThread =(endBlock-startBlock)%m_thread;
    for(size_t i = 0; i < maxThread; ++i){
        garbageCollection.push_back(new LinkageThread(correction, endBlock, &m_linkage, &m_linkageSqrt, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
    }

    for(size_t i = startBlock; i < endBlock; ++i){
        garbageCollection[i%maxThread]->Addstart(i); //So we distribute the items to the corresponding "thread"
    }
    std::vector<pthread_t*> threadList;
    for(size_t i = 0; i < maxThread; ++i){
        pthread_t *thread1 = new pthread_t();
        threadList.push_back(thread1);
        int threadStatus = pthread_create( thread1, NULL, &LinkageThread::triangularProcess, garbageCollection[i]);
        if(threadStatus != 0){
            throw "Failed to spawn thread with status: "+std::to_string(threadStatus);
        }
    }
    for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter){
		pthread_join(*threadList[threadIter], NULL);
    }

    for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter){
        delete(threadList[threadIter]);
	}
	threadList.clear();
    for(size_t cleaning = 0; cleaning < garbageCollection.size(); ++cleaning){
        delete(garbageCollection[cleaning]);
    }
    garbageCollection.clear();
}


void Linkage::rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype){
    std::vector<LinkageThread*> garbageCollection;
    std::vector<pthread_t*> threadList;
    size_t snpStart = start-height;
    size_t snpEnd = start;
    if(m_thread > (snpEnd-snpStart)){
        //if there are more threads than items, each will do one horizontal row
        for(size_t i = snpStart; i < snpEnd; ++i){
            garbageCollection.push_back(new LinkageThread(correction, i, i+1, start, width, &m_linkage, &m_linkageSqrt, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
            pthread_t *thread1 = new pthread_t();
            threadList.push_back(thread1);
            int threadStatus = pthread_create( thread1, NULL, &LinkageThread::rectangularProcess, garbageCollection.back());
            if(threadStatus != 0){
                throw "Failed to spawn thread with status: "+std::to_string(threadStatus);
            }
        }
    }
    else{
        size_t maxThread = (snpEnd -snpStart ) / m_thread;
        if(maxThread >=1) maxThread = m_thread;
        else maxThread =(snpEnd -snpStart ) % m_thread;
        size_t step=(snpEnd -snpStart )/maxThread;
        size_t remaining = (snpEnd -snpStart )%maxThread;
        size_t current = snpStart;
        for(size_t i = 0; i < maxThread; ++i){
            if(remaining > 0){
                garbageCollection.push_back(new LinkageThread(correction, current, current+step+1, start, width,&m_linkage, &m_linkageSqrt, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
                remaining--;
                current++;
            }
            else{
                garbageCollection.push_back(new LinkageThread(correction, current, current+step, start, width, &m_linkage, &m_linkageSqrt, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
            }
            pthread_t *thread1 = new pthread_t();
            threadList.push_back(thread1);
            int threadStatus = pthread_create( thread1, NULL, &LinkageThread::rectangularProcess, garbageCollection.back());
            if(threadStatus != 0){
                throw "Failed to spawn thread with status: "+std::to_string(threadStatus);
            }
            current += step;
        }
    }
    for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter){
		pthread_join(*threadList[threadIter], NULL);
    }

    for(size_t threadIter = 0; threadIter < threadList.size(); ++threadIter){
        delete(threadList[threadIter]);
	}
    for(size_t i = 0; i < garbageCollection.size(); ++i){
        delete garbageCollection[i];
    }
}

ProcessCode Linkage::Initialize(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize){
    m_perfectLd.clear();
	if(genotype.empty()){
        return continueProcess;
	}
	if(blockSize == 0){
        //Doesn't have to build the LD matrix when the block size is 0
        return continueProcess;
	}

	//Use previous informations
    if(prevResiduals == 0){
        m_linkage = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
        m_linkageSqrt = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
    }
    else{
		Eigen::MatrixXd temp = m_linkage.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
		temp = m_linkageSqrt.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkageSqrt= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkageSqrt.topLeftCorner(prevResiduals, prevResiduals) = temp;

    }
    return continueProcess;
}

void Linkage::Reinitialize(size_t &genotypeSize){
    m_perfectLd.clear();

	if(genotypeSize == 0){
        throw "This is not possible unless there are some complicated undetected bug. Please contact the author with the input\nNo genotype left after removing perfect LDs";

    }
    else if(genotypeSize < (unsigned) m_linkage.cols()){
        m_linkage.conservativeResize(genotypeSize, genotypeSize);
        m_linkageSqrt.conservativeResize(genotypeSize, genotypeSize);
        //m_varLinkage.conservativeResize(genotypeSize, genotypeSize);
    }
}

void Linkage::Construct(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize, bool correction){
    m_perfectLd.clear(); //this should have served its purpose. Clear and recycle.
	if(genotype.empty()){
        throw "There is no genotype to work on";
	}
	if(blockSize == 0){
        //Doesn't have to build the LD matrix when the block size is 0
        throw "Block size is 0, something must be wrong.";
	}

    size_t numOfSnp = genotype.size();
    //Construct the LD
	std::vector<int> startLoc;
	size_t currentBlockSize = blockSize;
	if(currentBlockSize > numOfSnp){ //The block contains all the available Snps. Should usually only happen when it is the last block.
        currentBlockSize = numOfSnp;
	}
    size_t stepSize = currentBlockSize/3; //We want to increment the step by third

    if(stepSize == 0){ //Only possible when there is less than 3 Snps (max, 2x2). Therefore don't bother using multi-thread
        for(size_t i = 0; i < numOfSnp; ++i){
            m_linkage(i,i) = 1.0;
            m_linkageSqrt(i,i) = 1.0;
            for(size_t j = i+1; j < numOfSnp; ++j){ //invert the direction so that we can do it faster
                double rSquare =0.0;
                double r = 0.0;
                //if(fabs((*m_snpList)[(*m_snpLoc)[i]]->Getbp()-(*m_snpList)[(*m_snpLoc)[j]]->Getbp()) <= 2000000){
                    genotype[i]->GetbothR(genotype[j], correction, r, rSquare);
                //}
                m_linkage(i,j) = rSquare;
                m_linkage(j,i) = rSquare;
                m_linkageSqrt(i,j) = r;
                m_linkageSqrt(j,i) = r;
            }
        }
    }
    else if(stepSize > 0){ //Trying to work using the multi-threading method
        size_t counting= 0;
        for(size_t i = 0; i < genotype.size(); i+= stepSize){
            if(counting < 2){
                if(prevResiduals==0){
                    startLoc.push_back(i);
                }
            }
            else{
                startLoc.push_back(i);
            }
            counting++;
        }
        if(!startLoc.empty() && genotype.size() - startLoc.back() < stepSize) startLoc.pop_back(); //If the remaining is less than a stepSize from the end
        if((startLoc.empty() || startLoc.size()==0) && genotype.size()!= 0 ){
            //When this is the only block left
            triangularThread(stepSize*2, genotype.size(), correction, genotype);
            rectangularThread(0, genotype.size(),stepSize*2, correction, genotype);
        }
        else{
            for(size_t i = 0; i  < startLoc.size()-1; ++i){
                triangularThread(startLoc[i], startLoc[i+1], correction, genotype);
            }
            if(startLoc.size() != 0 && !startLoc.empty()){ //Double safe
                triangularThread(startLoc.back(), genotype.size(), correction, genotype);
            }

            for(size_t i = 0; i < startLoc.size()-1; ++i){
                if(prevResiduals==0 && i == 0){

                } //do nothing
                else if(prevResiduals==0 && i==1){
                    //The smaller chunk
                    rectangularThread(startLoc[i], startLoc[i+1], stepSize, correction, genotype);
                }
                else{
                    rectangularThread(startLoc[i], startLoc[i+1], stepSize*2, correction, genotype);
                }
            }
            if(startLoc.size() >= 1){
                rectangularThread(startLoc.back(),genotype.size() , stepSize*2, correction, genotype);
            }
        }
        startLoc.clear();
    }
    else{
        throw "Undefined behaviour! Step size should never be negative as block size is positive";
	}

}


size_t Linkage::Remove(){
    if(m_perfectLd.empty()) return 0;
    else{
        //Do something stupid first to make it easier for me for now
        std::map<size_t, bool> requireRemove;
        for(size_t i=0; i < m_perfectLd.size(); ++i){
            requireRemove[m_perfectLd[i]] = true;
        }
        //Eigen::Matrix(row, col) (Just so I remember)
        size_t rowIndex = 0;
        size_t colIndex = 0;
        size_t numRow = m_linkage.rows();
        size_t numCol = m_linkage.cols();

        for(size_t i = 0; i < numRow; ++i){
            colIndex = rowIndex;
            if(requireRemove.find(i)==requireRemove.end()){
                for(size_t j = i; j < numCol; ++j){
                    if(requireRemove.find(j) == requireRemove.end()){
                        m_linkage(rowIndex, colIndex) = m_linkage(i, j);
                        m_linkage(colIndex,rowIndex ) = m_linkage(j, i);
                        m_linkageSqrt(rowIndex, colIndex) = m_linkageSqrt(i, j);
                        m_linkageSqrt(colIndex,rowIndex ) = m_linkageSqrt(j, i);
                        colIndex++;
                    }
                }
                for(size_t j = colIndex; j < numCol; ++j){
                    m_linkage(rowIndex, j) = 0.0;
                    m_linkage(j, rowIndex) = 0.0;
                    m_linkageSqrt(rowIndex, j) = 0.0;
                    m_linkageSqrt(j, rowIndex) = 0.0;
                }
                rowIndex++;
            }
        }
        for(size_t i = rowIndex; i < numRow; ++i){
            for(size_t j = i; j < numCol; ++j){
                m_linkage(i, j)  = 0.0;
                m_linkage(j, i) = 0.0;
                m_linkageSqrt(i, j)  = 0.0;
                m_linkageSqrt(j, i) = 0.0;
            }
        }
        return m_perfectLd.size();
    }
}

//Updating the two corresponding structures
void Linkage::Update(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc){
    size_t offset = 0;
    for(size_t i = 0; i < m_perfectLd.size(); ++i){
        Genotype *temp = genotype[m_perfectLd[i]-offset];
        genotype.erase(genotype.begin() + (m_perfectLd[i]-offset));
        delete temp;
        snpLoc.erase(snpLoc.begin() + (m_perfectLd[i]-offset));
        offset++;
    }

}

Eigen::VectorXd Linkage::solve(size_t start, size_t length, Eigen::VectorXd const *const betaEstimate, Eigen::VectorXd const *const sqrtChiSq, Eigen::VectorXd *perSnpEffect,  Eigen::VectorXd *effectiveReturnResult,  size_t sampleSize, size_t snpStart){
    /** Perform the eigen value decomposition here */
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m_linkage.block(start, start, length, length));
    /** Calculate the tolerance threshold */
    double tolerance = std::numeric_limits<double>::epsilon() * length * es.eigenvalues().array().maxCoeff();
    /** Generate the pseudo inverse by removing any eigenvalues less than the tolerance threshold */
    Eigen::MatrixXd rInverse = es.eigenvectors()*(es.eigenvalues().array() > tolerance).select(es.eigenvalues().array().inverse(), 0).matrix().asDiagonal() * es.eigenvectors().transpose();


    Eigen::VectorXd result= rInverse*(*betaEstimate).segment(start, length);
    /** Here we try to perform the iterative adjustment to reduce the relative error */
    Eigen::VectorXd error =m_linkage.block(start, start, length, length)*result - (*betaEstimate).segment(start, length);
	double bNorm = (*betaEstimate).segment(start, length).norm();
    double relative_error = error.norm() / bNorm;
    double prev_error = relative_error+1;
    Eigen::VectorXd update = result;
// TODO (swchoi#1#): Might actually want to add a hard termination here. E.g. terminate after X cycles ...
//

    while(relative_error < prev_error){
        prev_error = relative_error;
        update.noalias()=rInverse*(-error);
        relative_error = 0.0;
        error.noalias()= m_linkage.block(start, start, length, length)*(result+update) - (*betaEstimate).segment(start, length);
        relative_error = error.norm() / bNorm;
        if(relative_error < 1e-300) relative_error = 0;
        result = result+update;
    }

    Eigen::VectorXd effectiveNumber = Eigen::VectorXd::Constant(length,1.0);
    double eNorm = effectiveNumber.norm();
    Eigen::VectorXd effectiveResult = rInverse*effectiveNumber;
    error = m_linkage.block(start, start, length, length)*effectiveResult-effectiveNumber;
    relative_error = error.norm()/eNorm;
    prev_error = relative_error+1;
    update=effectiveResult;
    while(relative_error < prev_error){
        prev_error = relative_error;
        update.noalias()=rInverse*(-error);
        relative_error = 0.0;
        error.noalias()= m_linkage.block(start, start, length, length)*(effectiveResult+update) - effectiveNumber;
        relative_error = error.norm() / eNorm;
        if(relative_error < 1e-300) relative_error = 0;
        effectiveResult = effectiveResult+update;
    }

    (*effectiveReturnResult) = effectiveResult;
    (*perSnpEffect) =m_linkage.block(start, start, length, length).colwise().sum();
    /** Here we try to calculate the variance */
    /*
    Eigen::MatrixXd ncpEstimate = (4*m_linkageSqrt.block(start, start, length, length)).array()*((*sqrtChiSq).segment(start, length)*(*sqrtChiSq).segment(start, length).transpose()).array();
    Eigen::VectorXd minusF = (Eigen::VectorXd::Constant(length, 1.0)-(*betaEstimate).segment(start, length));
    for(size_t i = 0; i < length; ++i){
        minusF(i) = minusF(i)/((double)sampleSize-2.0+((*sqrtChiSq).segment(start, length))(i)*((*sqrtChiSq).segment(start, length))(i));
    }
    (*variance).noalias() = (rInverse*minusF.asDiagonal()*(ncpEstimate-2*m_linkage.block(start, start, length, length))*minusF.asDiagonal()*rInverse);
*/

    return result;
}



Eigen::MatrixXd Linkage::solve(size_t start, size_t length, Eigen::MatrixXd const *const sampleMatrix, size_t sampleSize, size_t snpStart){
    /** Perform the eigen value decomposition here */
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m_linkage.block(start, start, length, length));
    /** Calculate the tolerance threshold */
    double tolerance = std::numeric_limits<double>::epsilon() * length * es.eigenvalues().array().maxCoeff();
    /** Generate the pseudo inverse by removing any eigenvalues less than the tolerance threshold */
    Eigen::MatrixXd rInverse = es.eigenvectors()*(es.eigenvalues().array() > tolerance).select(es.eigenvalues().array().inverse(), 0).matrix().asDiagonal() * es.eigenvectors().transpose();


    Eigen::MatrixXd result= rInverse*(*sampleMatrix).block(start, 0,length, (*sampleMatrix).cols());
    /** Here we try to perform the iterative adjustment to reduce the relative error */
    Eigen::MatrixXd error =m_linkage.block(start, start, length, length)*result - (*sampleMatrix).block(start, 0,length, (*sampleMatrix).cols());
	double bNorm = (*sampleMatrix).block(start, 0,length, (*sampleMatrix).cols()).norm();
    double relative_error = error.norm() / bNorm;
    double prev_error = relative_error+1;
    Eigen::MatrixXd update = result;
// TODO (swchoi#1#): Might actually want to add a hard termination here. E.g. terminate after X cycles ...
//

    while(relative_error < prev_error){
        prev_error = relative_error;
        update.noalias()=rInverse*(-error);
        relative_error = 0.0;
        error.noalias()= m_linkage.block(start, start, length, length)*(result+update) - (*sampleMatrix).block(start, 0,length, (*sampleMatrix).cols());
        relative_error = error.norm() / bNorm;
        if(relative_error < 1e-300) relative_error = 0;
        result = result+update;
    }
    return result;
}


void Linkage::print(){
    std::ofstream testing;
    testing.open("TESTING");
    testing << m_linkageSqrt << std::endl;
    testing.close();
}
