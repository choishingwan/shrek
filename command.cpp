#include "command.h"

Command::Command(){}
Command::~Command(){}

void Command::initialize(int argc, char* argv[]){
    std::cerr << "------------------------------------------------------------------------------"  << std::endl;
    std::cerr << "| Snp HeRitability Estimation Kit                                             |" << std::endl;
    std::cerr << "| version "<<m_version <<"                                                                |" << std::endl;
    std::cerr << "| (C) 2014 Johnny Kwan, Sam Choi                                              |" << std::endl;
    std::cerr << "| The University of Hong Kong                                                 |" << std::endl;
    std::cerr << "| GNU General Public License v2.0                                             |" << std::endl;
    std::cerr << "------------------------------------------------------------------------------"  << std::endl;

    if(argc==1){
        std::cerr << "Usage: ./SHREK <command> [options]"                                      << std::endl;
        std::cerr << "Command:    quant          Quantitative Trait"                                   << std::endl;
        std::cerr << "            cc             Case control study"                                   << std::endl;
        std::cerr << "            risk-qt        Risk Prediction on continuous traits"                 << std::endl;
        std::cerr << "            risk-cc        Risk Prediction on dichotomize traits"                << std::endl;
        std::cerr << std::endl;
        std::cerr << "To see the specific parameters of each mode, use: "                              << std::endl;
        std::cerr << "       ./"<<argv[0]<<" <command> -h"                                             << std::endl;
        throw std::runtime_error("Unspecified mode");
    }
	m_programmeName =argv[0];
	std::string mode(argv[1]);
    if(mode.compare("quant")==0){
        m_qt=true;
        if(argc == 2) printQuantUsage();
        else quantitativeProcess(argc, argv);
    }
    else if(mode.compare("cc")==0){
        m_cc=true;
        if(argc==2) printCCUsage();
        else caseControlProcess(argc, argv);
    }
    else if(mode.compare("risk-qt")==0){
        m_rqt=true;
        if(argc==2) printRiskQtUsage();
        else continuousRiskProcess(argc, argv);
    }
    else if(mode.compare("risk_cc")==0){
        m_rcc = true;
        if(argc ==2) printRiskQtUsage();
        else dichotomusRiskProcess(argc, argv);
    }
    else{
        throw std::runtime_error("Unspecified mode");
    }
}


void Command::getIndex(const std::vector<std::string> &index, const std::string &pvalueFileName, std::vector<size_t> &indexResult){
    //check if we can open the file
    std::ifstream pfile;
    if(pvalueFileName.empty()){
        throw "P-value file name not provided!";
    }
    pfile.open(pvalueFileName.c_str());
    if(!pfile.is_open()){
        throw "Cannot open pvalue file";
    }
    std::string header;
    std::getline(pfile, header);
    pfile.close();
    header = usefulTools::trim(header);
    if(header.empty()){
        throw "Empty header line for pvalue file";
    }
    std::vector<std::string> token;
    usefulTools::tokenizer(header, "\t ", &token);
    if(token.size() < index.size()){
        throw "Not enough field for header, please note that we assume the delimiter as space or tab";
    }

    // Not very efficient here, but should be ok consider the number of possible header
    // Maybe someone who is interested in optimization can optimize this
    bool found = false;
    for(size_t i = 0; i < index.size(); ++i){
        for(size_t j = 0; j < token.size(); ++j){
            if(index[i].compare(token[j])==0){
                found = true;
                indexResult.push_back(j);
                break;
            }
        }
        if(!found){
            std::cerr << index[i] << " not found in header" << std::endl;
        }

        found = false;
    }

}

void Command::caseControlProcess(int argc, char* argv[]){

    static const char *optString = "a:R:k:p:b:c:r:l:s:o:f:t:D:d:L:vunh?";
    static const struct option longOpts[]={
	    //Parameters for correction
        {"case",required_argument, NULL,'a'},
        {"control",required_argument, NULL,'R'},
        {"prevalence",required_argument, NULL,'k'},
        //Parameters on file input
        {"pfile",required_argument, NULL,'p'},
        {"bfile",required_argument, NULL,'b'},
        {"chr",required_argument, NULL,'c'},
        {"rs",required_argument, NULL,'r'},
        {"loc",required_argument, NULL,'l'},
        {"stat",required_argument, NULL,'s'},
        {"direction",required_argument, NULL,'D'},
        //General parameters
        {"out",required_argument, NULL,'o'},
        {"maf",required_argument, NULL,'f'},
        {"thread",required_argument, NULL,'t'},
        {"distance",required_argument, NULL,'d'},
        {"region",required_argument, NULL,'L'},
        {"validate",no_argument, NULL,'v'},
        {"pvalue",no_argument, NULL,'u'},
        {"correct",no_argument, NULL,'n'},
		{"help",no_argument,NULL,'h'},
		{NULL, 0, 0, 0}
	};

	std::vector<std::string> colNamesForIndex;
	std::vector<char> typeForIndex;
	int longIndex=0;
	int opt = 0;
	opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    while(opt!=-1){
		switch(opt){
		    //CC specific parameters
			case 'a':
				m_caseSize= atoi(optarg);
				break;
			case 'R':
				m_controlSize= atoi(optarg);
				break;
			case 'k':
				m_prevalence = atof(optarg);
				m_providedPrevalence =true;
				break;
            //File parameters
			case 'p':
				m_pValueFileName = optarg;
				break;
			case 'b':
                m_ldFilePrefix = optarg;
                break;
			case 'c':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('c');
                //m_chrIndex = atoi(optarg)-1;
                break;
			case 'r':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('r');
				//m_rsIndex = atoi(optarg)-1;
				break;
			case 'l':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('b');
				//m_bpIndex=atoi(optarg)-1;
				break;
            case 's':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('s');
                //m_stats = atoi(optarg)-1;
                break;
            case 'D':
                colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('D');
        //General parameters
			case 'o':
				m_outputPrefix = optarg;
				break;
			case 'f':
                m_maf = atof(optarg);
                m_providedMaf= true;
                break;
			case 't':
				m_thread = atoi(optarg);
				break;
			case 'd':
				m_distance = atoi(optarg);
				break;
			case 'L':
				m_regionList = optarg;
                break;
			case 'v':
				m_validate = true;
				break;
			case 'u':
                m_isPvalue = true;
                break;
			case 'n':
				m_ldCorrection = false;
				break;
			case 'h':
			case '?':
 				printCCUsage();
                throw 0;
				break;
			default:
				throw "Undefined operator, please use --help for more information!";
		}
		opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    std::vector<size_t> index;
    getIndex(colNamesForIndex, m_pValueFileName, index);
    if(index.size() != colNamesForIndex.size()){
        throw "Some field(s) not found in header, please check";
    }
    for(size_t i = 0; i < typeForIndex.size(); ++i){
        switch(typeForIndex[i]){
            case 'c':
                m_chrIndex = index[i];
                break;
            case 'b':
                m_bpIndex = index[i];
                break;
            case 'r':
                m_rsIndex = index[i];
                break;
            case 's':
                m_stats = index[i];
                break;
            case 'D':
                m_dirIndex = index[i];
                m_dirGiven =true;
                break;
        }
    }
    //Check for duplication
    sort( index.begin(), index.end() );
    index.erase( unique( index.begin(), index.end() ), index.end() );
    sort( typeForIndex.begin(), typeForIndex.end() );
    typeForIndex.erase( unique( typeForIndex.begin(), typeForIndex.end() ), typeForIndex.end() );
    if(index.size() != typeForIndex.size()){
        throw "Duplicated fields, please check your input";
    }
}

void Command::quantitativeProcess(int argc, char* argv[]){
    static const char *optString = "e:N:x:b:p:c:r:l:s:d:D:f:L:o:t:nuvh?";
	static const struct option longOpts[]={
	    //Qt specific parameter
		{"extreme",required_argument,NULL,'e'},
        {"sampleSize",required_argument,NULL,'N'},
        {"sampleIndex",required_argument,NULL,'x'},
        //File parameters
		{"bfile",required_argument,NULL,'b'},
		{"pfile",required_argument,NULL,'p'},
		{"chr",required_argument,NULL,'c'},
		{"rs",required_argument,NULL,'r'},
		{"loc",required_argument,NULL,'l'},
		{"stats", required_argument, NULL, 's'},
		{"direction",required_argument,NULL,'D'},
        //General parameters
		{"distance",required_argument,NULL,'d'},
		{"maf",required_argument,NULL,'f'},
		{"region",required_argument,NULL,'L'},
		{"out",required_argument,NULL,'o'},
		{"thread",required_argument,NULL,'t'},
		{"correct",no_argument,NULL,'n'},
		{"pvalue",no_argument,NULL,'u'},
		{"validate",no_argument,NULL,'v'},
		{"help",no_argument,NULL,'h'},
		{NULL, 0, 0, 0}
	};

	std::vector<std::string> colNamesForIndex;
	std::vector<char> typeForIndex;
	int longIndex=0;
	int opt = 0;
	opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    while(opt!=-1){
		switch(opt){
		    //Qt specific parameter
			case 'e':
				m_extremeAdjust = atof(optarg);
				m_provideExtremeAdjustment = true;
				break;
			case 'N':
				m_sampleSize= atoi(optarg);
				m_provideSampleSize = true;
				break;
			case 'x':
				//m_sampleSizeIndex = atoi(optarg)-1;
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('a');
				break;
            //File parameters
			case 'b':
                m_ldFilePrefix = optarg;
                break;
			case 'p':
				m_pValueFileName = optarg;
				break;
			case 'c':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('c');
                //m_chrIndex = atoi(optarg)-1;
                break;
			case 'r':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('r');
				//m_rsIndex = atoi(optarg)-1;
				break;
			case 'l':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('b');
				//m_bpIndex=atoi(optarg)-1;
				break;
            case 'D':
                colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('D');
            case 's':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('s');
                //m_stats = atoi(optarg)-1;
                break;
        //General parameters
            case 'd':
                m_distance = atoi(optarg);
                break;
			case 'f':
                m_maf = atof(optarg);
                m_providedMaf= true;
                break;
			case 'L':
				m_regionList = optarg;
                break;
			case 'o':
				m_outputPrefix = optarg;
				break;
			case 't':
				m_thread = atoi(optarg);
				break;
			case 'n':
				m_ldCorrection = false;
				break;
			case 'u':
                m_isPvalue = true;
                break;
			case 'v':
				m_validate = true;
				break;
    		case 'h':
			case '?':
 				printQuantUsage();
                throw 0;
				break;
			default:
				throw "Undefined operator, please use --help for more information!";
		}
		opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    std::vector<size_t> index;
    getIndex(colNamesForIndex, m_pValueFileName, index);
    if(index.size() != colNamesForIndex.size()){
        throw "Some field(s) not found in header, please check";
    }
    for(size_t i = 0; i < typeForIndex.size(); ++i){
        switch(typeForIndex[i]){
            case 'c':
                m_chrIndex = index[i];
                break;
            case 'b':
                m_bpIndex = index[i];
                break;
            case 'r':
                m_rsIndex = index[i];
                break;
            case 's':
                m_stats = index[i];
                break;
            case 'a':
                m_sampleSizeIndex = index[i];
                break;
            case 'D':
                m_dirIndex = index[i];
                m_dirGiven = true;
                break;
        }
    }
    //Check for duplication
    sort( index.begin(), index.end() );
    index.erase( unique( index.begin(), index.end() ), index.end() );
    sort( typeForIndex.begin(), typeForIndex.end() );
    typeForIndex.erase( unique( typeForIndex.begin(), typeForIndex.end() ), typeForIndex.end() );
    if(index.size() != typeForIndex.size()){
        throw "Duplicated fields, please check your input";
    }
    if(!m_isPvalue && !m_dirGiven){ //Test statistic from Qt trait contains the sign
        m_dirGiven = true;
        m_dirIndex = m_stats;
    }
    //TODO Check input, e.g. sample size should all be integers, if input is p-value, then it should be in the range of 0 and 1
}

void Command::dichotomusRiskProcess(int argc, char* argv[]){
    static const char *optString = "E:T:D:g:Aa:R:k:p:b:c:r:l:so:f:t:d:L:vunh?";
	static const struct option longOpts[]={
        //Risk prediction specific parameters
        {"ref", required_argument, NULL, 'E'},
        {"alt", required_argument, NULL, 'T'},
        {"dir", required_argument, NULL, 'D'},
        {"geno", required_argument, NULL, 'g'},
        {"ambiguous", required_argument, NULL, 'A'},
	    //CC specific parameters
        {"case",required_argument, NULL,'a'},
        {"control",required_argument, NULL,'R'},
        {"prevalence",required_argument, NULL,'k'},
        //Parameters on file input
        {"pfile",required_argument, NULL,'p'},
        {"bfile",required_argument, NULL,'b'},
        {"chr",required_argument, NULL,'c'},
        {"rs",required_argument, NULL,'r'},
        {"loc",required_argument, NULL,'l'},
        {"stat",required_argument, NULL,'s'},
        //General parameters
        {"out",required_argument, NULL,'o'},
        {"maf",required_argument, NULL,'f'},
        {"thread",required_argument, NULL,'t'},
        {"distance",required_argument, NULL,'d'},
        {"region",required_argument, NULL,'L'},
        {"validate",no_argument, NULL,'v'},
        {"pvalue",no_argument, NULL,'u'},
        {"correct",no_argument, NULL,'n'},
		{"help",no_argument,NULL,'h'},
		{NULL, 0, 0, 0}
	};
	int longIndex=0;
	int opt = 0;
	opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    while(opt!=-1){
		switch(opt){
		    //risk prediction specific parameters
            case 'E':
                m_ref = atoi(optarg)-1;
                break;
            case 'T':
                m_alt = atoi(optarg)-1;
                break;
            case 'D':
                m_dirIndex = atoi(optarg)-1;
                break;
            case 'G':
                m_genotypeFilePrefix = optarg;
                break;
            case 'A':
                m_keep = true;
                break;

		    //CC specific parameters
			case 'a':
				m_caseSize= atoi(optarg);
				break;
			case 'R':
				m_controlSize= atoi(optarg);
				break;
			case 'k':
				m_prevalence = atof(optarg);
				m_providedPrevalence =true;
				break;
            //File parameters
			case 'p':
				m_pValueFileName = optarg;
				break;
			case 'b':
                m_ldFilePrefix = optarg;
                break;
			case 'c':
                m_chrIndex = atoi(optarg)-1;
                break;
			case 'r':
				m_rsIndex = atoi(optarg)-1;
				break;
			case 'l':
				m_bpIndex=atoi(optarg)-1;
				break;
            case 's':
                //m_stats = processRange(optarg);
                m_stats = atoi(optarg)-1;
                break;
        //General parameters
			case 'o':
				m_outputPrefix = optarg;
				break;
			case 'f':
                m_maf = atof(optarg);
                m_providedMaf= true;
                break;
			case 't':
				m_thread = atoi(optarg);
				break;
			case 'd':
				m_distance = atoi(optarg);
				break;
			case 'L':
				m_regionList = optarg;
                break;
			case 'v':
				m_validate = true;
				break;
			case 'u':
                m_isPvalue = true;
                break;
			case 'n':
				m_ldCorrection = false;
				break;
			case 'h':
			case '?':
 				printRiskCCUsage();
                throw 0;
				break;
			default:
				throw "Undefined operator, please use --help for more information!";
		}
		opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    }


}
void Command::continuousRiskProcess(int argc, char* argv[]){
    static const char *optString = "E:T:D:g:AN:x:b:p:c:r:l:s:d:f:L:o:t:nuvh?";
	static const struct option longOpts[]={
        //Risk prediction specific parameters
        {"ref", required_argument, NULL, 'E'},
        {"alt", required_argument, NULL, 'T'},
        {"dir", required_argument, NULL, 'D'},
        {"geno", required_argument, NULL, 'g'},
        {"ambiguous", required_argument, NULL, 'A'},

	    //Qt specific parameter
        {"sampleSize",required_argument,NULL,'N'},
        {"sampleIndex",required_argument,NULL,'x'},
        //File parameters
		{"bfile",required_argument,NULL,'b'},
		{"pfile",required_argument,NULL,'p'},
		{"chr",required_argument,NULL,'c'},
		{"rs",required_argument,NULL,'r'},
		{"loc",required_argument,NULL,'l'},
		{"stats", required_argument, NULL, 's'},
        //General parameters
		{"distance",required_argument,NULL,'d'},
		{"maf",required_argument,NULL,'f'},
		{"region",required_argument,NULL,'L'},
		{"out",required_argument,NULL,'o'},
		{"thread",required_argument,NULL,'t'},
		{"correct",no_argument,NULL,'n'},
		{"pvalue",no_argument,NULL,'u'},
		{"validate",no_argument,NULL,'v'},
		{"help",no_argument,NULL,'h'},
		{NULL, 0, 0, 0}
	};

	int longIndex=0;
	int opt = 0;
	opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    while(opt!=-1){
		switch(opt){
		    //risk prediction specific parameters
            case 'E':
                m_ref = atoi(optarg)-1;
                break;
            case 'T':
                m_alt = atoi(optarg)-1;
                break;
            case 'D':
                m_dirIndex = atoi(optarg)-1;
                break;
            case 'G':
                m_genotypeFilePrefix = optarg;
                break;
            case 'A':
                m_keep = true;
                break;

		    //Qt specific parameter
			case 'N':
				m_sampleSize= atoi(optarg);
				m_provideSampleSize = true;
				break;
			case 'x':
				m_sampleSizeIndex = atoi(optarg)-1;
				break;
            //File parameters
			case 'b':
                m_ldFilePrefix = optarg;
                break;
			case 'p':
				m_pValueFileName = optarg;
				break;
			case 'c':
                m_chrIndex = atoi(optarg)-1;
                break;
			case 'r':
				m_rsIndex = atoi(optarg)-1;
				break;
			case 'l':
				m_bpIndex=atoi(optarg)-1;
				break;
            case 's':
                //m_stats = processRange(optarg);
                m_stats = atoi(optarg)-1;
                break;
        //General parameters
            case 'd':
                m_distance = atoi(optarg);
                break;
			case 'f':
                m_maf = atof(optarg);
                m_providedMaf= true;
                break;
			case 'L':
				m_regionList = optarg;
                break;
			case 'o':
				m_outputPrefix = optarg;
				break;
			case 't':
				m_thread = atoi(optarg);
				break;
			case 'n':
				m_ldCorrection = false;
				break;
			case 'u':
                m_isPvalue = true;
                break;
			case 'v':
				m_validate = true;
				break;
    		case 'h':
			case '?':
 				printRiskQtUsage();
                throw 0;
				break;
			default:
				throw "Undefined operator, please use --help for more information!";
		}
		opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    }

}

std::vector<size_t> Command::processRange(std::string input){
    std::vector<size_t> index;
    if(input.empty()){
        throw std::runtime_error("Invalid index for statistic");
    }
    std::string temp = usefulTools::trim(input);
    std::vector<std::string> token;
    usefulTools::tokenizer(temp, ",", &token);
    for(size_t i =0; i < token.size(); ++i){
        if(usefulTools::trim(token[i]).empty()){
            throw std::runtime_error("Invalid index for statistic, no number between \",\"");
        }
        std::vector<std::string> range;
        usefulTools::tokenizer(token[i], "-", &range);

        if(range.size() == 1){
            index.push_back(atoi(range[0].c_str())-1);
        }
        else if(range.size() ==2){
            size_t start = atoi(range[0].c_str())-1;
            size_t end = atoi(range[1].c_str());
            for(size_t j = start; j < end; ++j){
                index.push_back(j);
            }
        }
        else{
            throw std::runtime_error("Invalid format for statistic index, should be <num>-<num>");
        }
    }
    //check if there are any duplication
    sort(index.begin(), index.end());
    size_t before = index.size();
    index.erase( unique( index.begin(), index.end() ), index.end());
    size_t after = index.size();
    if(before != after){
        throw std::runtime_error("Duplicated index for statistic");
    }
    return index;
}


void Command::printRunSummary(std::string regionMessage){
    std::cerr 	<< "SHREK\tCoding by Sam CHOI\tMethod by Johnny KWAN" <<std::endl
        << "===============================================================" << std::endl
		<< "Performing analysis using the following parameters: " << std::endl
		<< "===============================================================" << std::endl
		<< "Essential Input  " <<std::endl;
		std::cerr	<< "Linkage File Prefix  : " << m_ldFilePrefix << std::endl;
	    std::cerr 	<< "P-Value File         : " << m_pValueFileName << std::endl;
	    std::cerr   << "Input is P-Value     : ";
    if(m_isPvalue) std::cerr << "True" << std::endl;
    else std::cerr << "False" << std::endl;
    if(!m_outputPrefix.empty()){
        std::cerr   << "Output File          : " << m_outputPrefix << std::endl;
    }
    if(m_cc){
		std::cerr	<< "Mode                 : Case Control " << std::endl
					<< "Number of Case       : " << m_caseSize << std::endl
					<< "Number of Control    : " << m_controlSize << std::endl
					<< "Prevalence           : " << m_prevalence << std::endl << std::endl;
	}
	else if(m_qt){
		std::cerr	<< "Mode                 : Quantitative Trait" << std::endl;
		if(m_provideSampleSize) std::cerr	<< "Sample Size          : " << m_sampleSize << std::endl;
		if(m_provideExtremeAdjustment) std::cerr <<   "Extreme Adjustment   : " << m_extremeAdjust << std::endl << std::endl;
	}
	else if(m_rqt){
        std::cerr	<< "Mode                 : Risk Prediction (Continuous Trait)" << std::endl;
        std::cerr   << "Genotype File Prefix : " << m_genotypeFilePrefix << std::endl;
        std::cerr   << "Ref Allele           : " << m_ref << std::endl;
        std::cerr   << "Alt Allele           : " << m_alt << std::endl;
		if(m_provideSampleSize) std::cerr	<< "Sample Size          : " << m_sampleSize << std::endl;
	}
	else if(m_rcc){
        std::cerr	<< "Mode                 : Risk Prediction (Dichotomous Trait)" << std::endl;
        std::cerr   << "Genotype File Prefix : " << m_genotypeFilePrefix << std::endl;
        std::cerr   << "Ref Allele           : " << m_ref << std::endl;
        std::cerr   << "Alt Allele           : " << m_alt << std::endl
					<< "Number of Case       : " << m_caseSize << std::endl
					<< "Number of Control    : " << m_controlSize << std::endl
					<< "Prevalence           : " << m_prevalence << std::endl << std::endl;
		if(m_provideSampleSize) std::cerr	<< "Sample Size          : " << m_sampleSize << std::endl;
	}
    std::cerr	<< "===============================================================" << std::endl
				<< "Options " << std::endl
				<< "Number of Thread     : " << m_thread << std::endl;
    if(m_ldCorrection){
        std::cerr << "Use LD correction    : True" << std::endl;
    }
    else{
        std::cerr << "Use LD correction    : False" << std::endl;
    }
    std::cerr	<< "LD distance          : " << m_distance << std::endl;
	std::cerr	<< "Number of regions    : " << regionMessage << std::endl;

	std::cerr << std::endl << std::endl;

}



void Command::printUsage(){
    //General Usage
    std::cerr << "General options: " << std::endl;
    std::cerr << "  -d,--distance    The maximum distance allowed between SNPs at the start of "   << std::endl;
    std::cerr << "                   window and the end of the window. The larger the value, the " << std::endl;
    std::cerr << "                   slower the programme runs "                                   << std::endl;
    std::cerr << "  -f,--maf         The maf threshold for the reference genotype file. SNPs with "<< std::endl;
    std::cerr << "                   maf less than this value will not be used for the analysis"   << std::endl;
    std::cerr << "  -u,--pvalue      Indicate whether if the input is p-value or test-statistic "  << std::endl;
    std::cerr << "                   Cannot handle p-value of 0"                                   << std::endl;
    std::cerr << "  -n,--correct     Turn off LD correction. The LD correction is used when "      << std::endl;
    std::cerr << "                   using the genotype file to calculate the LD matrix. The "     << std::endl;
    std::cerr << "                   LD correction is performed to adjust for number of sample"    << std::endl;
    std::cerr << "                   used for calculating the LD. "                                << std::endl;
    std::cerr << "                   We uses Weir & Hill(1980)'s formula to correct for the bias. "<< std::endl;
    std::cerr << "  -L,--region      The region files. You may provide a bed file in the format:"  << std::endl;
    std::cerr << "                   <region name>:<fileName>,<region name>:<fileName>,..."        << std::endl;
    std::cerr << "                   The summary output will provide the per region estimate "     << std::endl;
    std::cerr << "                   where the region are label by their order of input"           << std::endl;
    std::cerr << "  -t,--thread      The number of thread use. This will affect not only the speed"<< std::endl;
    std::cerr << "                   but also the memory requirement of the programme as we will " << std::endl;
    std::cerr << "                   load more SNPs into the memory at one time if more threads "  << std::endl;
    std::cerr << "                   are available"                                                << std::endl;
    std::cerr << "  -v,--validate    Validate the SNPs. If this option is set, the programme will" << std::endl;
    std::cerr << "                   compare the SNPs information in the p-value file with those " << std::endl;
    std::cerr << "                   in the genotype file. Will output a warning if the "          << std::endl;
    std::cerr << "                   information does not match. This validation is still very "   << std::endl;
    std::cerr << "                   naive so it is important for the user to make sure the "      << std::endl;
    std::cerr << "                   information matched. "                                        << std::endl;
    std::cerr << "  -o,--out         The output file prefix. If provided, the programme will "     << std::endl;
    std::cerr << "                   generate a <output>.sum and <output>.res file where the "     << std::endl;
    std::cerr << "                   <output>.sum file will provide the run summary and the "      << std::endl;
    std::cerr << "                   <output>.res file will provide the detail statistics"         << std::endl;
    std::cerr << "  -h,-?,--help     Display the detail help message (This message)"               << std::endl;
    std::cerr << "                                                                               " << std::endl;
}
void Command::printCCUsage(){
    std::cerr << "Estimation of Heritability from Case Control Study:" << std::endl;
    std::cerr << "usage: ./shrek cc [options]"                                                     << std::endl;
    std::cerr << "Required options: "                                                              << std::endl;
    std::cerr << "  -p,--pfile       The p-value file. The test statistic for each snp must "      << std::endl;
    std::cerr << "                   be provided. "                                                << std::endl;
    std::cerr << "  -b,--bfile       The linkage  file prefix. The programme will use this to "    << std::endl;
    std::cerr << "                   calculate the LD matrix. Will require the fam, bim and bed"   << std::endl;
    std::cerr << "                   file. Please try to perform quality control beforehand "      << std::endl;
    std::cerr << "  -s,--stats       The column header of test statistic or p-value in the p-value"<< std::endl;
    std::cerr << "                   file. It must be provided for the programme to run"           << std::endl;
    std::cerr << "  -c,--chr         The column header of chromosome in the p-value file"          << std::endl;
    std::cerr << "  -r,--rs          The column header of rsid in the p-value file"                << std::endl;
    std::cerr << "  -l,--loc         The column header of rsid coordinate in the p-value file"     << std::endl;
    std::cerr << "  -k,--prevalence  Specify the prevalence of the phenotype"                      << std::endl;
    std::cerr << "  -a,--case        The number of case used in the study"                         << std::endl;
    std::cerr << "  -R,--control     The number of control used in the study"                      << std::endl;
    std::cerr << "                                                                               " << std::endl;
    printUsage();
    exit(0);
}
void Command::printQuantUsage(){

    std::cerr << "Estimation of Heritability from Quantitative Trait Study:" << std::endl;
    std::cerr << "usage: ./shrek quant [options]"                                                  << std::endl;
    std::cerr << "Required options: "                                                              << std::endl;
    std::cerr << "  -p,--pfile       The p-value file. The test statistic for each snp must "      << std::endl;
    std::cerr << "                   be provided. "                                                << std::endl;
    std::cerr << "  -b,--bfile       The linkage  file prefix. The programme will use this to "    << std::endl;
    std::cerr << "                   calculate the LD matrix. Will require the fam, bim and bed"   << std::endl;
    std::cerr << "                   file. Please try to perform quality control beforehand "      << std::endl;
    std::cerr << "  -s,--stats       The column header of test statistic or p-value in the p-value"<< std::endl;
    std::cerr << "                   file. It must be provided for the programme to run"           << std::endl;
    std::cerr << "  -c,--chr         The column header of chromosome in the p-value file"          << std::endl;
    std::cerr << "  -r,--rs          The column header of rsid in the p-value file"                << std::endl;
    std::cerr << "  -l,--loc         The column header of rsid coordinate in the p-value file"     << std::endl;
    std::cerr << "  -x,--sampleIndex The column header of sample size in the p-value file. "       << std::endl;
    std::cerr << "                   If not provided, one must provide the sample size using the"  << std::endl;
    std::cerr << "                   -N or --sampleSize option."                                   << std::endl;
    std::cerr << "  -N,--sampleSize  The number of sample used in the association study. Must be " << std::endl;
    std::cerr << "                   provided if the sampleIndex was not provided"                 << std::endl;
    std::cerr << "  -e,--extreme     The extreme phenotype adjustment value. Should be: "          << std::endl;
    std::cerr << "                   Variance after selection / Variance before selection"         << std::endl;
    std::cerr << "                                                                               " << std::endl;
    printUsage();
    exit(0);
}

void Command::printRiskCCUsage(){

}
void Command::printRiskQtUsage(){

}


