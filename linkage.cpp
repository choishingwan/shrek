#include "linkage.h"

<<<<<<< HEAD
<<<<<<< HEAD
Linkage::Linkage(std::vector<Snp*> *snpList, size_t thread):m_snpList(snpList), m_thread(thread){
	m_numItem= 0;
    m_effectiveNumber=0.0;

}
=======
Linkage::Linkage(size_t thread, std::vector<Snp*> *snpList, std::deque<size_t> *snpLoc):m_thread(thread), m_snpList(snpList), m_snpLoc(snpLoc){}
>>>>>>> perfectLd
=======
Linkage::Linkage(size_t thread, std::vector<Snp*> *snpList, std::deque<size_t> *snpLoc):m_thread(thread), m_snpList(snpList), m_snpLoc(snpLoc){}
>>>>>>> perfectLD

Linkage::~Linkage()
{
	//dtor
}

size_t Linkage::rows() const { return m_linkage.rows(); }
size_t Linkage::cols() const { return m_linkage.cols(); }
Eigen::MatrixXd Linkage::block(size_t blockStart, size_t lengthOfBlock){ return m_linkage.block(blockStart, blockStart, lengthOfBlock, lengthOfBlock); }

void Linkage::triangularThread( const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype, std::vector<size_t> &perfectLd, std::deque<size_t> &snpLoc){
    //Make the thread region
    std::vector<LinkageThread*> garbageCollection;
    size_t maxThread = (endBlock-startBlock)/ m_thread;
    if(maxThread >=1) maxThread = m_thread;
    else maxThread =(endBlock-startBlock)%m_thread;
    for(size_t i = 0; i < maxThread; ++i){
<<<<<<< HEAD
<<<<<<< HEAD
        garbageCollection.push_back(new LinkageThread(correction, endBlock, &m_linkage, &genotype, &snpLoc, m_snpList, &perfectLd));
=======
        garbageCollection.push_back(new LinkageThread(correction, endBlock, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
>>>>>>> perfectLd
=======
        garbageCollection.push_back(new LinkageThread(correction, endBlock, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
>>>>>>> perfectLD
    }

    for(size_t i = startBlock; i < endBlock; ++i){
        garbageCollection[i%maxThread]->Addstart(i); //So we distribute the items to the corresponding "thread"
    }
    std::vector<pthread_t*> threadList;
    for(size_t i = 0; i < maxThread; ++i){
        pthread_t *thread1 = new pthread_t();
        threadList.push_back(thread1);
        std::cerr << "Create thread" << std::endl;
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


void Linkage::rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype, std::vector<size_t> &perfectLd, std::deque<size_t> &snpLoc){
    std::vector<LinkageThread*> garbageCollection;
    std::vector<pthread_t*> threadList;
    size_t snpStart = start-height;
    size_t snpEnd = start;
    if(m_thread > (snpEnd-snpStart)){
        //if there are more threads than items, each will do one horizontal row
        for(size_t i = snpStart; i < snpEnd; ++i){
<<<<<<< HEAD
<<<<<<< HEAD
            garbageCollection.push_back(new LinkageThread(correction, i, i+1, start, width, &m_linkage, &genotype, &snpLoc, m_snpList, &perfectLd));
=======
            garbageCollection.push_back(new LinkageThread(correction, i, i+1, start, width, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
>>>>>>> perfectLd
=======
            garbageCollection.push_back(new LinkageThread(correction, i, i+1, start, width, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
>>>>>>> perfectLD
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
<<<<<<< HEAD
<<<<<<< HEAD
                garbageCollection.push_back(new LinkageThread(correction, current, current+step+1, start, width, &m_linkage, &genotype, &snpLoc, m_snpList, &perfectLd));
=======
                garbageCollection.push_back(new LinkageThread(correction, current, current+step+1, start, width, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
>>>>>>> perfectLd
=======
                garbageCollection.push_back(new LinkageThread(correction, current, current+step+1, start, width, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
>>>>>>> perfectLD
                remaining--;
                current++;
            }
            else{
<<<<<<< HEAD
<<<<<<< HEAD
                garbageCollection.push_back(new LinkageThread(correction, current, current+step, start, width, &m_linkage, &genotype, &snpLoc, m_snpList, &perfectLd));
=======
                garbageCollection.push_back(new LinkageThread(correction, current, current+step, start, width, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
>>>>>>> perfectLd
=======
                garbageCollection.push_back(new LinkageThread(correction, current, current+step, start, width, &m_linkage, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
>>>>>>> perfectLD
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

<<<<<<< HEAD
<<<<<<< HEAD

void Linkage::Remove(std::vector<size_t> &perfectLd, std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc){
    //First remove by rows
    size_t offset = 0;
    if(perfectLd.empty()) return;
    std::cerr << "we start the removal process" << std::endl;
    for(size_t i = 0; i < perfectLd.size(); ++i){
            std::cerr << "Remove " << i << " row and column" << std::endl;
        size_t rowToRemove = perfectLd[i]-offset;
        size_t numRows = m_linkage.rows()-1;
        size_t numCols = m_linkage.cols();
        m_linkage.block(rowToRemove,0,numRows-rowToRemove,numCols) = m_linkage.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
        m_linkage.conservativeResize(numRows,numCols);
        numRows = m_linkage.rows();
        numCols = m_linkage.cols()-1;
        size_t colToRemove = perfectLd[i]-offset;
        m_linkage.block(0,colToRemove,numRows,numCols-colToRemove) = m_linkage.block(0,colToRemove+1,numRows,numCols-colToRemove);
        m_linkage.conservativeResize(numRows,numCols);
        offset++;
    }
    offset = 0;
    for(size_t i = 0; i < perfectLd.size(); ++i){
        genotype.erase(genotype.begin() + (perfectLd[i]-offset));
        snpLoc.erase(snpLoc.begin() + (snpLoc[i]-offset));
        offset++;
    }
}

ProcessCode Linkage::Initialize(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize){
=======
ProcessCode Linkage::Initialize(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize){
    m_perfectLd.clear();
>>>>>>> perfectLd
=======
ProcessCode Linkage::Initialize(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize){
    m_perfectLd.clear();
>>>>>>> perfectLD
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
    else if(prevResiduals < (unsigned) m_linkage.cols()){ //Reasonable size
		Eigen::MatrixXd temp = m_linkage.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
<<<<<<< HEAD
<<<<<<< HEAD
        return continueProcess;
    }
    else{
        std::cerr << "Should not happen. Definitely something went wrong with the algorithm. Please contact the author if you see this message" << std::endl;
        std::cerr << "Error message: prevResiduals larger than m_linkage size" << std::endl;
        return fatalError;
=======
=======
>>>>>>> perfectLD
    }
    return continueProcess;
}

ProcessCode Linkage::Reinitialize(size_t &genotypeSize){
    m_perfectLd.clear();

	if(genotypeSize == 0){
        std::cerr << "This is not possible unless there are some complicated undetected bug. Please contact the author with the input" << std::endl;
        std::cerr << "No genotype left after removing perfect LDs" << std::endl;
        return fatalError;
<<<<<<< HEAD
    }
    else if(genotypeSize < (unsigned) m_linkage.cols()){
        m_linkage.conservativeResize(genotypeSize, genotypeSize);
>>>>>>> perfectLd
    }
    return continueProcess;
}

<<<<<<< HEAD
ProcessCode Linkage::Reinitialize(const size_t genotypeSize){
	if(genotypeSize==0){
        return continueProcess;
	}
    m_linkage.conservativeResize(genotypeSize, genotypeSize);
    return continueProcess;
}

ProcessCode Linkage::Construct(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize, bool correction, std::vector<size_t> &perfectLd, std::deque<size_t> &snpLoc){
=======
ProcessCode Linkage::Construct(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize, bool correction){
    m_perfectLd.clear();
>>>>>>> perfectLd
=======
    }
    else if(genotypeSize < (unsigned) m_linkage.cols()){
        m_linkage.conservativeResize(genotypeSize, genotypeSize);
    }
    return continueProcess;
}

ProcessCode Linkage::Construct(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize, bool correction){
    m_perfectLd.clear();
>>>>>>> perfectLD
	if(genotype.empty()){
        return continueProcess;
	}
	if(blockSize == 0){
        //Doesn't have to build the LD matrix when the block size is 0
        return continueProcess;
	}

<<<<<<< HEAD
<<<<<<< HEAD
=======

>>>>>>> perfectLd
=======

>>>>>>> perfectLD
    //Construct the LD
	std::vector<size_t> startLoc;
	size_t currentBlockSize = blockSize;
	if(currentBlockSize > genotype.size()){
        currentBlockSize = genotype.size();
	}
    size_t stepSize = currentBlockSize/3;
    std::cerr << "Step size is: " << stepSize << std::endl;
    if(stepSize == 0){ //Doesn't require threading
            std::cerr << "No threading" << std::endl;
        for(size_t i = 0; i < genotype.size(); ++i){
            m_linkage(i,i) = 1.0;
<<<<<<< HEAD
<<<<<<< HEAD
            for(size_t j = genotype.size()-1 ; j >=i+1; --j){
                if(m_linkage(i,j) == 0.0){
                    double rSquare = genotype[i]->Getr(genotype[j], correction);
                    if(rSquare>=1.0){
                        perfectLd.push_back(j); //It is possible that the perfectLD isn't unique
                        (*m_snpList)[snpLoc[j]]->shareHeritability((*m_snpList)[snpLoc[i]]);
=======
=======
>>>>>>> perfectLD
            for(size_t j = genotype.size()-1; j > i; --j){ //invert the direction
                if(m_linkage(i,j) == 0.0){
                    double rSquare = genotype[i]->Getr(genotype[j], correction);
                    if(i != j && std::fabs(rSquare-1.0) < G_EPSILON_DBL){
                        m_perfectLd.push_back(j);
                        (*m_snpList)[(*m_snpLoc)[j]]->shareHeritability((*m_snpList)[(*m_snpLoc)[i]]);
<<<<<<< HEAD
>>>>>>> perfectLd
=======
>>>>>>> perfectLD
                    }
                    m_linkage(i,j) = rSquare;
                    m_linkage(j,i) = rSquare;
                }
<<<<<<< HEAD
<<<<<<< HEAD
                else break;
=======
                else break; //Already calculated before
>>>>>>> perfectLd
=======
                else break; //Already calculated before
>>>>>>> perfectLD
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
        std::cerr << "Perform threading: " << std::endl;
        if((startLoc.empty() || startLoc.size()==0) && genotype.size()!= 0 ){
            //When this is the only block left
            std::cerr << "Last block" << std::endl;
            triangularThread(stepSize*2, genotype.size(), correction, genotype, perfectLd, snpLoc);
            rectangularThread(0, genotype.size(),stepSize*2, correction, genotype, perfectLd, snpLoc);
        }
        else{
                std::cerr << "Before triangle" << std::endl;
            for(size_t i = 0; i  < startLoc.size()-1; ++i){
                triangularThread(startLoc[i], startLoc[i+1], correction, genotype, perfectLd, snpLoc);
            }
            if(startLoc.size() != 0 && !startLoc.empty()){ //Double safe
                triangularThread(startLoc.back(), genotype.size(), correction, genotype, perfectLd, snpLoc);
            }
            std::cerr << "After triangle" << std::endl;
            for(size_t i = 0; i < startLoc.size()-1; ++i){
                if(prevResiduals==0 && i == 0){

                } //do nothing
                else if(prevResiduals==0 && i==1){
                    //The smaller chunk
                    rectangularThread(startLoc[i], startLoc[i+1], stepSize, correction, genotype, perfectLd, snpLoc);
                }
                else{
                    rectangularThread(startLoc[i], startLoc[i+1], stepSize*2, correction, genotype, perfectLd, snpLoc);
                }
            }
            if(startLoc.size() >= 1){
                rectangularThread(startLoc.back(),genotype.size() , stepSize*2, correction, genotype, perfectLd, snpLoc);
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


Eigen::VectorXd Linkage::solve(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Eigen::VectorXd *effective){
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(m_linkage.block(start, start, length, length), Eigen::ComputeThinU);

    Eigen::MatrixXd rInvert = svd.matrixU()*(svd.singularValues().array().abs() > svd.threshold()).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().transpose();
    Eigen::VectorXd result = rInvert*(*betaEstimate).segment(start, length);
	double relative_error = 0.0;
	Eigen::VectorXd error =m_linkage.block(start, start, length, length)*result - (*betaEstimate).segment(start, length);
    relative_error = error.norm() / (*betaEstimate).segment(start, length).norm();

    double prev_error = relative_error+1;
    Eigen::VectorXd update = result;
    while(relative_error < prev_error){
        prev_error = relative_error;
        update =rInvert*(-(error));
        relative_error = 0.0;
        error= m_linkage.block(start, start, length, length)*(result+update) - (*betaEstimate).segment(start, length);
        relative_error = error.norm() / (*betaEstimate).segment(start, length).norm();
        //if(relative_error < 1e-300) relative_error = 0;
        result = result+update;
    }

    Eigen::VectorXd ones = Eigen::VectorXd::Constant(length, 1.0);
    (*effective) = rInvert*ones;
    error =m_linkage.block(start, start, length, length)*(*effective) - ones;
    relative_error = error.norm()/ones.norm();
    prev_error = relative_error+1;
    while(relative_error < prev_error){
        prev_error = relative_error;
        update =rInvert*(-(error));
        relative_error = 0.0;
        error= m_linkage.block(start, start, length, length)*((*effective)+update) - ones;
        relative_error = error.norm() / ones.norm();
        //if(relative_error < 1e-300) relative_error = 0;
        (*effective) = (*effective)+update;
    }
    return result;
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
<<<<<<< HEAD

}

<<<<<<< HEAD
void Linkage::print(){ //DEBUG
=======

void Linkage::print(){
>>>>>>> perfectLd
=======

}


void Linkage::print(){
>>>>>>> perfectLD
    std::cout << m_linkage << std::endl;
}
