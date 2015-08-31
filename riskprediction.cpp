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
    targetGenotype = new GenotypeFileHandler(commander->GetgenotypeFile(), commander->Getthread(), commander->GetoutputPrefix());
}


RiskPrediction::~RiskPrediction()
{
    //dtor
    delete targetGenotype;
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
                m_samplePheno.push_back(0.0);
            }
        }
    }
    famFile.close();
    size_t prevLoc = 0; //This is use to check if the ordering is correct. If the prev > current, then the ordering of the two file is different, and will cause problem.
    std::string prevChr = "";
    std::map<std::string, bool> dupCheck;
    duplicate = 0;
    size_t ambiguous=0;
    size_t problem = 0;
    size_t corDiff = 0;
    size_t positive = 0;
    size_t notFound = 0;
    std::map<std::string, size_t> concordanceIndex;
    size_t lineNum = 0;
    while(std::getline(bimFile, line)){
        line = usefulTools::trim(line);
        std::vector<std::string> token;
        usefulTools::tokenizer(line, "\t ", &token);
        if(!line.empty() && token.size() >= 6 ){
                lineNum++;
            std::string chr = token[0];
            std::string rsId = token[1];
            size_t loc = atoi(token[3].c_str());
            std::string refAllele = token[4];
            std::string altAllele = token[5];
            m_genoInclude.push_back(-1);
            m_flipCheck.push_back(false);
            if(snpIndex.find(rsId)!= snpIndex.end()){
                if(prevChr.empty()){
                    prevChr = chr;
                }
                else if(prevChr.compare(chr) != 0){
                    prevChr = chr;
                }
                else if(prevLoc > (*m_snpList)[snpIndex[rsId]]->Getbp()){
                    //std::cerr << prevLoc << "\t" << (*m_snpList)[snpIndex[rsId]]->Getbp() << std::endl;
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
                        int compareStatus = compareAllele(refAllele, altAllele, curRef, curAlt);
                        if(compareStatus==1){
                            //ok
                            if(isAmbiguous && !m_keep){
                                ambiguous++;
                            }
                            else{
                                if(isAmbiguous) ambiguous++;
                                m_flipCheck.back()=false;
                                dupCheck[rsId] = true;
                                positive++;
                                m_genoInclude.back() = sIndex;
                                concordanceIndex[rsId] = sIndex;
                            }
                        }
                        else if(compareStatus==2){
                            //flip
                            if(isAmbiguous && !m_keep){
                                ambiguous++;
                            }
                            else{
                                if(isAmbiguous) ambiguous++;
                                m_flipCheck.back()=true;
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
    //size_t totalNum = genotypeFileHandler->GetestimateSnpTotal()*3;
	while(process != completed && process != fatalError){
        //Now everything should almost be the same as that in the snpestimation except for the function called
        process = genotypeFileHandler->getSnps(genotype, snpLoc, m_snpList, chromosomeStart, chromosomeEnd, m_maf,prevResidual, blockSize);

		if(process == completed && !chromosomeEnd){

		}
		else{
			//Now calculate the LD matrix
			linkageMatrix->Initialize(genotype, prevResidual, blockSize);
			linkageMatrix->Construct(genotype, prevResidual, blockSize, m_ldCorrection);
            targetGenotype->Getsamples(&normalizedGenotype, snpLoc, m_snpList, genotype.size() - normalizedGenotype.rows(), m_flipCheck);
            numProcessed+= genotype.size(); //Finished the LD construction
            //Now perform the decomposition, which only involve distributing the matrix to different stuff to work on
            decompositionHandler->Decompose(blockSize, snpLoc, genotype, chromosomeStart, chromosomeEnd, m_samplePheno, normalizedGenotype);


            if(!chromosomeEnd){
				if(blockSize > genotype.size()) throw "When block size is bigger than the number of genotype, it must be the end of chromosome";
				size_t retain = blockSize/3*2;
				//Will also need to remove the front of the normalizedGenotype
                Eigen::MatrixXd temp = normalizedGenotype.bottomRows(blockSize/3*2);
                normalizedGenotype = temp;
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

void RiskPrediction::result(){
    if(!m_outPrefix.empty()){
        std::string outputName = m_outPrefix;
        outputName.append(".res");
        std::ofstream result;
        result.open(outputName.c_str());
        if(!result.is_open()){
            std::cerr << "Cannot open file: " << outputName << std::endl;
            std::cerr << "Will output to stdout" << std::endl;
            std::cout << "sampleId\tpheno"<<std::endl;
            for(size_t i = 0; i < m_sampleId.size(); ++i){
                std::cout << m_sampleId[i] << "\t" << m_samplePheno[i] << std::endl;
            }
        }
        else{
            result <<"sampleId\tpheno" << std::endl;
            for(size_t i = 0; i < m_sampleId.size(); ++i){
                result << m_sampleId[i] << "\t" << m_samplePheno[i] << std::endl;
            }
        }
    }
    else{
        std::cout << "sampleId\tpheno"<<std::endl;
        for(size_t i = 0; i < m_sampleId.size(); ++i){
            std::cout << m_sampleId[i] << "\t" << m_samplePheno[i] << std::endl;
        }
    }
}

std::string RiskPrediction::convert(std::string allele){
    if(allele.compare("A")==0 ||allele.compare("a")==0  ) return "T";
    if(allele.compare("T")==0 ||allele.compare("t")==0  ) return "A";
    if(allele.compare("C")==0 ||allele.compare("c")==0  ) return "G";
    if(allele.compare("G")==0 ||allele.compare("g")==0  ) return "C";
    return "-";
}

std::string RiskPrediction::upper(std::string &str){
    std::string converted;
	for(size_t i = 0; i < str.size(); ++i)
		converted += std::toupper(str[i]);

	return converted;
}
int RiskPrediction::compareAllele(std::string refAllele, std::string altAllele, std::string targetRef, std::string targetAlt){
    //same = 1
    //flip = 2;
    //different = 3;
    std::string refA = upper(refAllele);
    std::string altA = upper(altAllele);
    std::string tRefA = upper(targetRef);
    std::string tAltA = upper(targetAlt);
    if(refA.compare(tRefA)==0 && altA.compare(tAltA)==0){
        return 1;
    }
    else if(refA.compare(tAltA)==0 && altA.compare(tRefA)==0){
        return 2;
    }
    else if(convert(refA).compare(tRefA)==0 && convert(altA).compare(tAltA)==0){
        return 1;
    }
    else if(convert(refA).compare(tAltA)==0 && convert(altA).compare(tRefA)==0){
        return 2;
    }
    else return 3;
}

