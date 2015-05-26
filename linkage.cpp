#include "linkage.h"

std::mutex Linkage::mtx; //DEBUG
Linkage::Linkage(){
    m_perfectLd =std::vector<size_t>();
    m_snpLoc = nullptr;
    m_snpList = nullptr;
    m_thread = 1;
}
/*
Linkage::Linkage(size_t thread, std::vector<Snp*> *snpList, std::deque<size_t> *snpLoc):m_thread(thread), m_snpList(snpList), m_snpLoc(snpLoc){
    m_perfectLd=std::vector<size_t>();
}
*/
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

double Linkage::ExpectedR2(double rSq, size_t numSample, size_t predictor){
    return 0;
    /*gsl_sf_result result;
    int status = gsl_sf_hyperg_2F1_e(1, 1, 0.5 * (numSample + 1), rSq, &result);
    if(status == GSL_SUCCESS) {
        double y = result.val;
        double value = 1.0 - ((numSample - predictor- 1.0)/(numSample - 1.0)) * (1.0 - rSq) * y;
        (value>0) ? value=value : value=0.0;
        (value<1) ? value=value : value=1.0;
        return value;
    }
    else{
        std::cerr << "Problem with calculating the variance of Rsq" << std::endl;
        exit(-1);
    }
    */
}

double Linkage::VarianceR2(double rSq, size_t numSample, size_t predictor){
    return 0; //Deactivate the variance of R2 at the moment to avoid problem
    /*
    gsl_sf_result result;
    int status = gsl_sf_hyperg_2F1_e(2, 2, 0.5 * (numSample + 3), rSq, &result);
    if(status == GSL_SUCCESS){
        double y = result.val;
        double expected =(ExpectedR2(rSq,numSample,predictor) - 1);
        double value =(((numSample -predictor -1.0)*((double)numSample-predictor+1.0))/((double)numSample*(double)numSample-1.0)) * ((1.0 - rSq)*(1.0 - rSq)) *y-(expected*expected);
        return value;
    }
    else{
        std::cerr << "Problem with calculating the variance of Rsq" << std::endl;
        exit(-1);
    }
    */
}

size_t Linkage::rows() const { return m_linkage.rows(); }
size_t Linkage::cols() const { return m_linkage.cols(); }
Eigen::MatrixXd Linkage::block(size_t blockStart, size_t lengthOfBlock){ return m_linkage.block(blockStart, blockStart, lengthOfBlock, lengthOfBlock); }
Eigen::MatrixXd Linkage::varBlock(size_t blockStart, size_t lengthOfBlock){ return m_varLinkage.block(blockStart, blockStart, lengthOfBlock, lengthOfBlock); }
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
            garbageCollection.push_back(new LinkageThread(correction, i, i+1, start, width, &m_linkage, &m_linkageSqrt, &genotype, m_snpLoc, &m_perfectLd, m_snpList));
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
        m_linkageSqrt = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
        //m_varLinkage = Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
    }
    else{
		Eigen::MatrixXd temp = m_linkage.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
		temp = m_linkageSqrt.bottomRightCorner(prevResiduals,prevResiduals);
		m_linkageSqrt= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		m_linkageSqrt.topLeftCorner(prevResiduals, prevResiduals) = temp;

		//temp = m_varLinkage.bottomRightCorner(prevResiduals,prevResiduals);
		//m_varLinkage= Eigen::MatrixXd::Zero(genotype.size(), genotype.size());
		//m_varLinkage.topLeftCorner(prevResiduals, prevResiduals) = temp;
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
        m_linkageSqrt.conservativeResize(genotypeSize, genotypeSize);
        //m_varLinkage.conservativeResize(genotypeSize, genotypeSize);
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
            m_linkage(i,i) = 1.0; //When linkage is self, there is no variance
            m_linkageSqrt(i,i) = 1.0;
            //m_varLinkage(i,i) = Linkage::VarianceR2(1.0,genotype[i]->GetnumSample(),0);
            for(size_t j = genotype.size()-1; j > i; --j){ //invert the direction
                if(m_linkage(i,j) == 0.0){
                    double rSquare = genotype[i]->GetrSq(genotype[j], correction);
                    double r = genotype[i]->Getr(genotype[j], correction);
                    if(i != j && std::fabs(rSquare-1.0) < G_EPSILON_DBL){
                        m_perfectLd.push_back(j);
                        (*m_snpList)[(*m_snpLoc)[j]]->shareHeritability((*m_snpList)[(*m_snpLoc)[i]]);
                    }
                    //m_varLinkage(i,j) = Linkage::VarianceR2(rSquare,numSample,1);
                    //m_varLinkage(j,i) = m_varLinkage(i,j);
                    m_linkage(i,j) = rSquare;
                    m_linkage(j,i) = rSquare;
                    m_linkageSqrt(i,j) = r;
                    m_linkageSqrt(j,i) = r;
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
                        m_linkageSqrt(rowIndex, colIndex) = m_linkageSqrt(i, j);
                        m_linkageSqrt(colIndex,rowIndex ) = m_linkageSqrt(j, i);
                        //m_varLinkage(rowIndex, colIndex) = m_varLinkage(i, j);
                        //m_varLinkage(colIndex,rowIndex ) = m_varLinkage(j, i);
                        colIndex++;
                    }
                }
                for(size_t j = colIndex; j < numCol; ++j){
                    m_linkage(rowIndex, j) = 0.0;
                    m_linkage(j, rowIndex) = 0.0;
                    m_linkageSqrt(rowIndex, j) = 0.0;
                    m_linkageSqrt(j, rowIndex) = 0.0;
                    //m_varLinkage(rowIndex, j) = 0.0;
                    //m_varLinkage(j, rowIndex) = 0.0;
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
                //m_varLinkage(i, j)  = 0.0;
                //m_varLinkage(j, i) = 0.0;
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


Eigen::VectorXd Linkage::solveChi(size_t start, size_t length, Eigen::VectorXd const *const betaEstimate, Eigen::MatrixXd &variance, size_t sampleSize){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m_linkage.block(start, start, length, length));
    double tolerance = std::numeric_limits<double>::epsilon() * length * es.eigenvalues().array().maxCoeff();
    Eigen::MatrixXd rInverse = es.eigenvectors()*(es.eigenvalues().array() > tolerance).select(es.eigenvalues().array().inverse(), 0).matrix().asDiagonal() * es.eigenvectors().transpose();
    Eigen::VectorXd result= rInverse*(*betaEstimate).segment(start, length);
    Eigen::VectorXd error =m_linkage.block(start, start, length, length)*result - (*betaEstimate).segment(start, length);



	double bNorm = (*betaEstimate).segment(start, length).norm();
    double relative_error = error.norm() / bNorm;
    double prev_error = relative_error+1;
    Eigen::VectorXd update = result;
    while(relative_error < prev_error){
        prev_error = relative_error;
        update=rInverse*(-error);
        relative_error = 0.0;
        error= m_linkage.block(start, start, length, length)*(result+update) - (*betaEstimate).segment(start, length);
        relative_error = error.norm() / bNorm;
        if(relative_error < 1e-300) relative_error = 0;
        result = result+update;
    }

    Eigen::VectorXd minusF = Eigen::VectorXd::Constant(length, 1.0)-(*betaEstimate).segment(start, length);
    Eigen::VectorXcd complexF =Eigen::VectorXcd::Zero(length);
    for(size_t i =0; i < length; ++i){
        complexF(i) = sqrt(((*betaEstimate).segment(start, length))(i));
    }
    variance = (rInverse*(((sampleSize-2)/((sampleSize*sampleSize-1)*(sampleSize-1)))*(minusF.asDiagonal()*(2*m_linkage.block(start, start, length, length).cast<std::complex<double> >()+4*(complexF.asDiagonal()*m_linkageSqrt.block(start, start, length, length).cast<std::complex<double> >()*complexF.adjoint().asDiagonal())*((sampleSize-2)*(sampleSize-1)+4)/(sampleSize+3))*minusF.asDiagonal()))*rInverse).real();

    return result;
}

Eigen::VectorXd Linkage::solve(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Eigen::VectorXd *effective){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m_linkage.block(start, start, length, length));
    double tolerance = std::numeric_limits<double>::epsilon() * length * es.eigenvalues().array().abs().maxCoeff();
    Eigen::MatrixXd rInverse = (es.eigenvalues().array() > tolerance).select(es.eigenvalues().array().sqrt().inverse(), 0).matrix().asDiagonal() * es.eigenvectors().transpose();
    Eigen::MatrixXd r = es.eigenvectors()*(es.eigenvalues().array() > tolerance).select(es.eigenvalues().array().sqrt(), 0).matrix().asDiagonal();
    Eigen::VectorXd result= rInverse*(*betaEstimate).segment(start, length);
    double relative_error = 0.0;
	Eigen::VectorXd error =r*result - (*betaEstimate).segment(start, length);
    relative_error = error.norm() / (*betaEstimate).segment(start, length).norm();
    double prev_error = relative_error+1;
    Eigen::VectorXd update = result;
    while(relative_error < prev_error){
        prev_error = relative_error;

        update=rInverse*(-error);
        relative_error = 0.0;
        error= r*(result+update) - (*betaEstimate).segment(start, length);
        relative_error = error.norm() / (*betaEstimate).segment(start, length).norm();
        if(relative_error < 1e-300) relative_error = 0;
        result = result+update;
    }

    Eigen::VectorXd ones = Eigen::VectorXd::Constant(length, 1.0);
    (*effective)=rInverse*ones;
    error =r*(*effective) - ones;
    relative_error = error.norm()/ones.norm();
    prev_error = relative_error+1;
    while(relative_error < prev_error){
        prev_error = relative_error;
        update=rInverse*(-error);
        relative_error = 0.0;
        error=r*((*effective)+update) - ones;
        relative_error = error.norm() / ones.norm();
        if(relative_error < 1e-300) relative_error = 0;
        (*effective) = (*effective)+update;
    }
    return result;
}

