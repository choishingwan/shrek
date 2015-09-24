#include "command.h"

Command::Command(){}

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
        quantitativeProcess(argc, argv);
    }
    else if(mode.compare("cc")==0){
        m_cc=true;
        caseControlProcess(argc, argv);
    }
    else if(mode.compare("risk-qt")==0){
        m_rqt=true;
        continuousRiskProcess(argc, argv);
    }
    else if(mode.compare("risk_cc")==0){
        m_rcc = true;
        dichotomusRiskProcess(argc, argv);
    }
}

void Command::caseControlProcess(int argc, char* argv[]){

    static const char *optString = "a:R:k:p:b:c:r:l:so:f:t:M:m:d:L:vunh?";
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
        //General parameters
        {"out",required_argument, NULL,'o'},
        {"maf",required_argument, NULL,'f'},
        {"thread",required_argument, NULL,'t'},
        {"maxBlock",required_argument, NULL,'M'},
        {"minBlock",required_argument, NULL,'m'},
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
			case 'M':
				m_maxBlock = atoi(optarg);
				m_maxBlockSet =true;
				break;
			case 'm':
				m_minBlock = atoi(optarg);
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
				m_ldCorrection = true;
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


}
void Command::quantitativeProcess(int argc, char* argv[]){
    static const char *optString = "e:N:x:b:p:c:r:l:s:d:f:L:M:m:o:t:nuvh?";
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
        //General parameters
		{"distance",required_argument,NULL,'d'},
		{"maf",required_argument,NULL,'f'},
		{"region",required_argument,NULL,'L'},
		{"maxBlock",required_argument,NULL,'M'},
		{"minBlock",required_argument,NULL,'m'},
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
                m_distance = atoi(optarg)-1;
                break;
			case 'f':
                m_maf = atof(optarg);
                m_providedMaf= true;
                break;
			case 'L':
				m_regionList = optarg;
                break;
			case 'M':
				m_maxBlock = atoi(optarg);
				m_maxBlockSet =true;
				break;
			case 'm':
				m_minBlock = atoi(optarg);
				break;
			case 'o':
				m_outputPrefix = optarg;
				break;
			case 't':
				m_thread = atoi(optarg);
				break;
			case 'n':
				m_ldCorrection = true;
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


}
void Command::dichotomusRiskProcess(int argc, char* argv[]){
    static const char *optString = "E:T:D:g:Aa:R:k:p:b:c:r:l:so:f:t:M:m:d:L:vunh?";
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
        {"maxBlock",required_argument, NULL,'M'},
        {"minBlock",required_argument, NULL,'m'},
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
			case 'M':
				m_maxBlock = atoi(optarg);
				m_maxBlockSet =true;
				break;
			case 'm':
				m_minBlock = atoi(optarg);
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
				m_ldCorrection = true;
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


}
void Command::continuousRiskProcess(int argc, char* argv[]){
    static const char *optString = "E:T:D:g:AN:x:b:p:c:r:l:s:d:f:L:M:m:o:t:nuvh?";
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
		{"maxBlock",required_argument,NULL,'M'},
		{"minBlock",required_argument,NULL,'m'},
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
                m_distance = atoi(optarg)-1;
                break;
			case 'f':
                m_maf = atof(optarg);
                m_providedMaf= true;
                break;
			case 'L':
				m_regionList = optarg;
                break;
			case 'M':
				m_maxBlock = atoi(optarg);
				m_maxBlockSet =true;
				break;
			case 'm':
				m_minBlock = atoi(optarg);
				break;
			case 'o':
				m_outputPrefix = optarg;
				break;
			case 't':
				m_thread = atoi(optarg);
				break;
			case 'n':
				m_ldCorrection = true;
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
        std::cerr   << "Ref Allele           : " << m_refIndex << std::endl;
        std::cerr   << "Alt Allele           : " << m_altIndex << std::endl;
		if(m_provideSampleSize) std::cerr	<< "Sample Size          : " << m_sampleSize << std::endl;
	}
	else if(m_rcc){
        std::cerr	<< "Mode                 : Risk Prediction (Dichotomous Trait)" << std::endl;
        std::cerr   << "Genotype File Prefix : " << m_genotypeFilePrefix << std::endl;
        std::cerr   << "Ref Allele           : " << m_refIndex << std::endl;
        std::cerr   << "Alt Allele           : " << m_altIndex << std::endl
					<< "Number of Case       : " << m_caseSize << std::endl
					<< "Number of Control    : " << m_controlSize << std::endl
					<< "Prevalence           : " << m_prevalence << std::endl << std::endl;
		if(m_provideSampleSize) std::cerr	<< "Sample Size          : " << m_sampleSize << std::endl;
	}
    std::cerr	<< "===============================================================" << std::endl
				<< "Options " << std::endl
				<< "Number of Thread     : " << m_thread << std::endl;
	if(m_maxBlock != 0){
        std::cerr << "Maximum block size   : " << m_maxBlock << std::endl;
	}
    if(m_minBlock != 0){
        std::cerr << "Minimum block size   : " << m_minBlock << std::endl;
    }
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



void printUsage(){

}
void printCCUsage(){

}
void printQuantUsage(){

}
void printRiskUsage(){

}


