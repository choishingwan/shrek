#include "linkage.h"


Linkage::Linkage(size_t thread, std::vector<Snp*> *snpList, std::deque<size_t> *snpLoc):m_thread(thread), m_snpList(snpList), m_snpLoc(snpLoc){}

Linkage::~Linkage()
{
	//dtor
}

size_t Linkage::rows() const { return m_linkage.rows(); }
size_t Linkage::cols() const { return m_linkage.cols(); }
Eigen::MatrixXd Linkage::block(size_t blockStart, size_t lengthOfBlock){ return m_linkage.block(blockStart, blockStart, lengthOfBlock, lengthOfBlock); }

void Linkage::triangularThread( const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype){
    //Make the thread region
    std::vector<LinkageThread*> garbageCollection;
    size_t maxThread = (endBlock-startBlock)/ m_thread;
    if(maxThread >=1) maxThread = m_thread;
    else maxThread =(endBlock-startBlock)%m_thread;
    for(size_t i = 0; i < maxThread; ++i){
        garbageCollection.push_back(new LinkageThread(correction, endBlock, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
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
            std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
            exit(-1);
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
            garbageCollection.push_back(new LinkageThread(correction, i, i+1, start, width, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
            pthread_t *thread1 = new pthread_t();
            threadList.push_back(thread1);
            int threadStatus = pthread_create( thread1, NULL, &LinkageThread::rectangularProcess, garbageCollection.back());
            if(threadStatus != 0){
                std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
                exit(-1);
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
                garbageCollection.push_back(new LinkageThread(correction, current, current+step+1, start, width, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
                remaining--;
                current++;
            }
            else{
                garbageCollection.push_back(new LinkageThread(correction, current, current+step, start, width, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
            }
            pthread_t *thread1 = new pthread_t();
            threadList.push_back(thread1);
            int threadStatus = pthread_create( thread1, NULL, &LinkageThread::rectangularProcess, garbageCollection.back());
            if(threadStatus != 0){
                std::cerr << "Failed to spawn thread with status: " << threadStatus << std::endl;
                exit(-1);
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
    }
    else{
		Eigen::MatrixXd temp = m_linkage.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
    }
    return continueProcess;
}

ProcessCode Linkage::Reinitialize(size_t &genotypeSize){
    m_perfectLd.clear();

	if(genotypeSize == 0){
        std::cerr << "This is not possible unless there are some complicated undetected bug. Please contact the author with the input" << std::endl;
        std::cerr << "No genotype left after removing perfect LDs" << std::endl;
        return fatalError;
    }
    else if(genotypeSize < (unsigned) m_linkage.cols()){
        m_linkage.conservativeResize(genotypeSize, genotypeSize);
    }
    return continueProcess;
}

ProcessCode Linkage::Construct(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize, bool correction){
    m_perfectLd.clear();
	if(genotype.empty()){
        return continueProcess;
	}
	if(blockSize == 0){
        //Doesn't have to build the LD matrix when the block size is 0
        return continueProcess;
	}


    //Construct the LD
	std::vector<int> startLoc;
	size_t currentBlockSize = blockSize;
	if(currentBlockSize > genotype.size()){
        currentBlockSize = genotype.size();
	}
    size_t stepSize = currentBlockSize/3;
    //If stepsize == 0 -> Block size < 3, just do the whole thing directly

    if(stepSize == 0){
        for(size_t i = 0; i < genotype.size(); ++i){
            m_linkage(i,i) = 1.0;
            for(size_t j = genotype.size()-1; j > i; --j){ //invert the direction
                if(m_linkage(i,j) == 0.0){
                    double rSquare = genotype[i]->Getr(genotype[j], correction);
                    if(i != j && std::fabs(rSquare-1.0) < G_EPSILON_DBL){
                        m_perfectLd.push_back(j);
                        (*m_snpList)[(*m_snpLoc)[j]]->shareHeritability((*m_snpList)[(*m_snpLoc)[i]]);
                    }
                    m_linkage(i,j) = rSquare;
                    m_linkage(j,i) = rSquare;
                }
                else break; //Already calculated before
            }
        }
    }
    else if(stepSize > 0){
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
        std::cerr << "Undefined behaviour! Step size should never be negative as block size is positive" << std::endl;
        return fatalError;
	}

    std::sort(m_perfectLd.begin(), m_perfectLd.end());
    m_perfectLd.erase( std::unique( m_perfectLd.begin(), m_perfectLd.end() ), m_perfectLd.end() );
	return completed;
}

size_t Linkage::Remove(){
    if(m_perfectLd.empty()) return 0;
    else{
            //Do something stupid first to make it easier for me for now
        std::map<size_t, bool> requireRemove;
        for(size_t i=0; i < m_perfectLd.size(); ++i){
            requireRemove[m_perfectLd[i]] = true;
        }
        //Eigen::Matrix(row, col)
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
                        colIndex++;
                    }
                }
                for(size_t j = colIndex; j < numCol; ++j){
                    m_linkage(rowIndex, j) = 0.0;
                    m_linkage(j, rowIndex) = 0.0;
                }
                rowIndex++;
            }
        }
        for(size_t i = rowIndex; i < numRow; ++i){
            for(size_t j = i; j < numCol; ++j){
                m_linkage(i, j)  = 0.0;
                m_linkage(j, i) = 0.0;
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

void Linkage::print(){ //DEBUG
    std::cout << m_linkage << std::endl;
}



Eigen::VectorXd Linkage::solve(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Eigen::VectorXd *effective){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m_linkage.block(start, start, length, length));
    double tolerance = std::numeric_limits<double>::epsilon() * length * es.eigenvalues().array().abs().maxCoeff();
    Eigen::MatrixXd A =es.eigenvectors()*(es.eigenvalues().array().abs() > tolerance).select(es.eigenvalues().array(), 0).matrix().asDiagonal() * es.eigenvectors().transpose(); //This should give us a well conditioned matrix
    Eigen::LDLT<Eigen::MatrixXd> ldlt(A); //Cholesky
    tolerance = std::numeric_limits<double>::epsilon() * length * ldlt.vectorD().array().abs().maxCoeff(); //This is to avoid having sqrt of negative number (or very small number);
    Eigen::MatrixXd D = (ldlt.vectorD().array()>tolerance).select(ldlt.vectorD().array().sqrt(), 0).matrix().asDiagonal();
    Eigen::MatrixXd ll = (ldlt.matrixL()*D);


    Eigen::VectorXd result=(*betaEstimate).segment(start, length);
    ll.triangularView<Eigen::Lower>().solveInPlace(result);
    if(result.maxCoeff() > 0.56){
        std::cout << m_linkage.block(start, start, length, length) << std::endl;
        std::cerr <<(*betaEstimate).segment(start, length) << std::endl;
        exit(-1);

    }


    return result;
    double relative_error = 0.0;
    //Eigen::VectorXd betaBar = es.eigenvectors()*es.eigenvectors().transpose()*(*betaEstimate).segment(start, length);
    Eigen::VectorXd betaBar =(*betaEstimate).segment(start, length); //Keep this simple first
    Eigen::VectorXd error =ll*result - betaBar;
    relative_error = error.norm() / betaBar.norm();

    double prev_error = relative_error+1;
    Eigen::VectorXd update = result;
    while(relative_error < prev_error){
        prev_error = relative_error;
        update = -error;
        ll.triangularView<Eigen::Lower>().solveInPlace(update);
        relative_error = 0.0;
        error= ll*(result+update) - betaBar;
        relative_error = error.norm() / betaBar.norm();
        if(relative_error < 1e-300) relative_error = 0;
        result = result+update;
    }

    Eigen::VectorXd ones = Eigen::VectorXd::Constant(length, 1.0);
    (*effective) = ones;
    ll.triangularView<Eigen::Lower>().solveInPlace((*effective));
    error =ll*(*effective) - ones;
    relative_error = error.norm()/ones.norm();
    prev_error = relative_error+1;
    while(relative_error < prev_error){
        prev_error = relative_error;
        update =(-(error));
        ll.triangularView<Eigen::Lower>().solveInPlace(update);
        relative_error = 0.0;
        error= ll*((*effective)+update) - ones;
        relative_error = error.norm() / ones.norm();
        if(relative_error < 1e-300) relative_error = 0;
        (*effective) = (*effective)+update;
    }
    return result;
}

Eigen::VectorXd Linkage::quickSolve(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Eigen::VectorXd *effective){
    Eigen::LDLT<Eigen::MatrixXd> ldlt(m_linkage.block(start, start, length, length));
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m_linkage.block(start, start, length, length));
    //Eigen::MatrixXd check = ldlt.matrixL();
    //std::cout << check << std::endl;
    //exit(-1);
    double tolerance = std::numeric_limits<double>::epsilon() * length * ldlt.vectorD().array().abs().maxCoeff();
    Eigen::MatrixXd D = (ldlt.vectorD().array()>tolerance).select(ldlt.vectorD().array().sqrt(), 0).matrix().asDiagonal();
    Eigen::MatrixXd ll = (ldlt.matrixL()*D).transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ll);
    Eigen::VectorXd result=Eigen::VectorXd::Zero(length);
    if(es.info() != Eigen::ComputationInfo::Success){
        std::cerr << "Failed to compute eigendecomposition" << std::endl;
    }
    else{
        //Build the eigenvalue
        tolerance = std::numeric_limits<double>::epsilon() * length * es.eigenvalues().array().abs().maxCoeff();
        Eigen::MatrixXd rInvert = es.eigenvectors()*(es.eigenvalues().array().abs() > tolerance).select(es.eigenvalues().array().inverse(), 0).matrix().asDiagonal() * es.eigenvectors().transpose();
        result = rInvert*(*betaEstimate).segment(start, length);
        double relative_error = 0.0;
        Eigen::VectorXd betaBar = es.eigenvectors()*es.eigenvectors().transpose()*(*betaEstimate).segment(start, length);
        //Eigen::VectorXd error =m_linkage.block(start, start, length, length)*result - (*betaEstimate).segment(start, length);
        //Eigen::VectorXd error =m_linkage.block(start, start, length, length)*result - betaBar;
        Eigen::VectorXd error =ll*result - betaBar;
        //relative_error = error.norm() / (*betaEstimate).segment(start, length).norm();
        relative_error = error.norm() / betaBar.norm();

        double prev_error = relative_error+1;
        Eigen::VectorXd update = result;
        while(relative_error < prev_error){
            prev_error = relative_error;
            update =rInvert*(-(error));
            relative_error = 0.0;
            //error= m_linkage.block(start, start, length, length)*(result+update) - (*betaEstimate).segment(start, length);
            //error= m_linkage.block(start, start, length, length)*(result+update) - betaBar;
            error= ll*(result+update) - betaBar;
            //relative_error = error.norm() / (*betaEstimate).segment(start, length).norm();
            relative_error = error.norm() / betaBar.norm();
            if(relative_error < 1e-300) relative_error = 0;
            result = result+update;
        }
        //Eigen::VectorXd ones = Eigen::VectorXd::Constant(length, 1.0);
        Eigen::VectorXd ones = es.eigenvectors()*es.eigenvectors().transpose()*Eigen::VectorXd::Constant(length, 1.0);
        (*effective) = rInvert*ones;
        //error =m_linkage.block(start, start, length, length)*(*effective) - ones;
        error =ll*(*effective) - ones;
        relative_error = error.norm()/ones.norm();
        prev_error = relative_error+1;
        while(relative_error < prev_error){
            prev_error = relative_error;
            update =rInvert*(-(error));
            relative_error = 0.0;
            //error= m_linkage.block(start, start, length, length)*((*effective)+update) - ones;
            error= ll*((*effective)+update) - ones;
            relative_error = error.norm() / ones.norm();
            if(relative_error < 1e-300) relative_error = 0;
            (*effective) = (*effective)+update;
        }
    }
    return result;
}
