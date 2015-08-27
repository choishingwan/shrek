#include "riskprediction.h"


RiskPrediction::RiskPrediction(const Command *commander,std::vector<Snp*> *snpList):m_snpList(snpList){
    m_thread = commander->Getthread();
	m_minBlock = commander->GetminBlock();
    m_maxBlock  = commander->GetmaxBlock();
    m_maf = commander->Getmaf();
    m_ldCorrection = commander->ldCorrect();
    m_keep = commander->keep();
    m_maxBlockSet = commander->maxBlockSet();
    m_genotypeFilePrefix = commander->GetgenotypeFile();
    m_ldFilePrefix = commander->GetldFilePrefix();
    m_outPrefix = commander->GetoutputPrefix();
    m_validate = commander->validate();
    targetGenotype = new GenotypeFileHandler(m_genotypeFilePrefix, m_thread, m_outPrefix);
}


RiskPrediction::~RiskPrediction()
{
    //dtor
}


void RiskPrediction::checkGenotype(){
    //First, get the SNP index
	size_t duplicate = 0;
	for(size_t i = 0; i < (*m_snpList).size(); ++i){
        if(snpIndex.find((*m_snpList)[i]->GetrsId())==snpIndex.end()){
            snpIndex.insert(std::pair<std::string,size_t>((*m_snpList)[i]->GetrsId(),i));
            (*m_snpList)[i]->addFlag(false);
        }
        else{
            duplicate++;
        }

    }
    if(duplicate == 0) std::cerr << "There are no duplicated rsID in the p-value file" << std::endl << std::endl;
    else std::cerr <<  "There are a total of " << duplicate << " duplicated rsID(s) in the p-value file" << std::endl << std::endl;

    std::string famFileName =m_genotypeFilePrefix;
    famFileName.append(".fam");
    std::ifstream famFile;
    famFile.open(famFileName.c_str());
    if(!famFile.is_open()){
        std::string message = "Cannot open fam file: ";
        message.append(famFileName);
        throw message;
    }

    std::string bimFileName = m_genotypeFilePrefix;
    bimFileName.append(".bim");
    std::ifstream bimFile;
    bimFile.open(bimFileName.c_str());
    if(!bimFile.is_open()){
        std::string message = "Cannot open bim file: ";
        message.append(bimFileName);
        throw message;
    }
    std::string line;
    while(std::getline(famFile, line)){
        line =usefulTools::trim(line);
        if(!line.empty()){
            std::vector<std::string> token;
            usefulTools::tokenizer(line, "\t ", &token);
            if(token.size() > 2){ //This will not be a proper fam file format, but as we only need the sample id, this is all what we need
                m_sampleId.push_back(token[1]);
            }
        }
    }
    famFile.close();
    size_t prevLoc = 0; //This is use to check if the ordering is correct. If the prev > current, then the ordering of the two file is different, and will cause problem.
    std::map<std::string, bool> dupCheck;
    duplicate = 0;
    size_t ambiguous=0;
    size_t problem = 0;
    size_t corDiff = 0;
    size_t positive = 0;
    size_t notFound = 0;
    std::map<std::string, size_t> concordanceIndex;
    while(std::getline(bimFile, line)){
        line = usefulTools::trim(line);
        std::vector<std::string> token;
        usefulTools::tokenizer(line, "\t ", &token);
        if(!line.empty() && token.size() >= 6 ){
            std::string chr = token[0];
            std::string rsId = token[1];
            size_t loc = atoi(token[3].c_str());
            std::string refAllele = token[4];
            std::string altAllele = token[5];
            m_genoInclude.push_back(-1);
            if(snpIndex.find(rsId)!= snpIndex.end()){
                if(prevLoc > (*m_snpList)[snpIndex[rsId]]->Getbp()){
                    throw "Ordering of the genotype file and the p-value file differ, please check if both file are coordinately sorted";
                }
                prevLoc =(*m_snpList)[snpIndex[rsId]]->Getbp();
                if(dupCheck.find(rsId)== dupCheck.end()){

                    size_t sIndex = snpIndex[rsId];
                    bool concordant = (*m_snpList)[sIndex]->Concordant(chr, loc, rsId);
                    if(!concordant){
                        corDiff++;
                    }
                    else{
                        //Now check the allele
                        std::string curRef = (*m_snpList)[sIndex]->Getref();
                        std::string curAlt = (*m_snpList)[sIndex]->Getalt();
                        bool isAmbiguous = Snp::ambiguousAllele(refAllele, altAllele);
                        //std::cerr << curRef << "\t" << curAlt << "\t" << refAllele << "\t" << altAllele << std::endl;
                        if(curRef.compare(refAllele)==0 && curAlt.compare(altAllele)==0){
                            //ok
                            if(isAmbiguous && !m_keep){
                                ambiguous++;
                            }
                            else{
                                if(isAmbiguous) ambiguous++;
                                m_flipCheck.push_back(false);
                                dupCheck[rsId] = true;
                                positive++;
                                m_genoInclude.back() = sIndex;
                                concordanceIndex[rsId] = sIndex;
                            }
                        }
                        else if(curRef.compare(altAllele)==0 && curAlt.compare(refAllele)==0){
                            //flip
                            if(isAmbiguous && !m_keep){
                                ambiguous++;
                            }
                            else{
                                if(isAmbiguous) ambiguous++;
                                m_flipCheck.push_back(true);
                                dupCheck[rsId] = true;
                                positive++;
                                m_genoInclude.back() = sIndex;
                                concordanceIndex[rsId] = sIndex;
                            }
                        }
                        else{
                            problem++;
                        }
                    }
                }
                else duplicate++;
            }
            else    notFound ++;
        }
    }
    bimFile.close();
    targetGenotype->initialize();
    size_t mafFilter = targetGenotype->mafCheck(m_genoInclude, m_sampleId.size());
    if(notFound != 0)    std::cerr << notFound << " SNPs does not have statistic information." << std::endl;
    if(duplicate != 0)    std::cerr << duplicate << " SNPs were duplicated in the genotype file." << std::endl;
    if(corDiff != 0)    std::cerr << corDiff << " SNPs with different coordinate removed." << std::endl;
    if(problem != 0)    std::cerr << problem << " SNPs with unresolved allele information removed." << std::endl;
    if(!m_keep && ambiguous != 0) std::cerr << ambiguous <<  " SNPs with ambiguous allele information removed." << std::endl;
    else if(ambiguous != 0) std::cerr << ambiguous <<  " SNPs with ambiguous allele information kept." << std::endl;
    if(mafFilter != 0) std::cerr << mafFilter <<  " monomorphic SNPs removed." << std::endl;
    std::cerr << positive  << " SNPs will be used for risk prediction" << std::endl;
    snpIndex = concordanceIndex;
}

void RiskPrediction::run(){
    GenotypeFileHandler *genotypeFileHandler = new GenotypeFileHandler(m_ldFilePrefix, m_thread, m_outPrefix);
    genotypeFileHandler->initialize(snpIndex, m_snpList, m_validate, m_maxBlockSet, m_maxBlock, m_minBlock,m_maf);
    targetGenotype->initialize();
    //To know whether if the SNP is in all three file, we only have to check the flag 0 of each SNP.
    //If false, then it is not included in the LD file.
	Genotype::SetsampleNum(genotypeFileHandler->GetsampleSize());
	ProcessCode process = startProcess;
    std::deque<Genotype*> genotype;
    Eigen::MatrixXd normalizedGenotype;
    std::deque<size_t> snpLoc;
    bool chromosomeStart= true;
    bool chromosomeEnd = false;
    size_t prevResidual;
    size_t blockSize;
    Linkage *linkageMatrix = new Linkage();
    linkageMatrix->setSnpList(m_snpList);
    linkageMatrix->setSnpLoc(&snpLoc);
    linkageMatrix->setThread(m_thread);
    Decomposition *decompositionHandler = new Decomposition( m_snpList, linkageMatrix, m_thread);
    size_t numProcessed = 0;
    size_t totalNum = genotypeFileHandler->GetestimateSnpTotal()*3;
	while(process != completed && process != fatalError){
        //Now everything should almost be the same as that in the snpestimation except for the function called
        process = genotypeFileHandler->getSnps(genotype, snpLoc, m_snpList, chromosomeStart, chromosomeEnd, m_maf,prevResidual, blockSize);

		if(process == completed && !chromosomeEnd){

		}
		else{
			//Now calculate the LD matrix
			linkageMatrix->Initialize(genotype, prevResidual, blockSize);
			linkageMatrix->Construct(genotype, prevResidual, blockSize, m_ldCorrection);
            targetGenotype->Getsamples(&normalizedGenotype, snpLoc, m_snpList, genotype.size() - normalizedGenotype.rows());

            numProcessed+= genotype.size(); //Finished the LD construction
            if(!chromosomeEnd){
				if(blockSize > genotype.size()) throw "When block size is bigger than the number of genotype, it must be the end of chromosome";
				size_t retain = blockSize/3*2;
				//Will also need to remove the front of the normalizedGenotype
                normalizedGenotype = normalizedGenotype.bottomRows(blockSize/3*2);
				Genotype::clean(genotype, retain);
				size_t removeCount = snpLoc.size() - retain;
				for(size_t i = 0; i < removeCount; ++i)	snpLoc.pop_front();

				numProcessed-=retain;
            }
            else{
                Genotype::clean(genotype,0);
                normalizedGenotype.resize(0,0);
                snpLoc.clear();
            }
        }
        if(chromosomeStart && !chromosomeEnd){
            chromosomeStart =false;
        }
        else if(chromosomeEnd){
            chromosomeStart = true;
            chromosomeEnd = false;
        }
        //Now we need to get the corresponding genotype based on snpLoc
        //This is the new thing, get the XB matrix
        //normalizedGenotype should contain previous stuff.

	}


    delete genotypeFileHandler;
    delete linkageMatrix;
}





