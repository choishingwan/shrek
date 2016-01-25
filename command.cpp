#include "command.h"

Command::Command(){}
Command::~Command(){}

void Command::initialize(int argc, char* argv[]){
    std::cerr << std::endl;
    std::cerr << "Program: SHREK (Snp HeRitability Estimation Kit)" << std::endl;
    std::cerr << "Version: "<< m_version << std::endl<<std::endl;
    // First check if we have provided any mode information, if not, then we will provide the general
    // error message
	std::string mode(argv[1]);
	// The concept here is that all the data checking and stuff are done within each individual functions
    // Including the error check
    if(mode.compare("quant")==0){
        m_qt=true;
        quantitativeProcess(argc, argv);
    }
    else if(mode.compare("binary")==0){
        m_cc=true;
        caseControlProcess(argc, argv);
    }
    /*
    else if(mode.compare("risk-qt")==0){
        m_rqt=true;
        continuousRiskProcess(argc, argv);
    }
    else if(mode.compare("risk-bin")==0){
        m_rcc = true;
        dichotomusRiskProcess(argc, argv);
    }
    */
    else if(mode.compare("help")==0 || mode.compare("--help")==0){
        //Provide the mode help message (redundant)
        std::cerr << "Usage: ./SHREK <command> [options]"                                         << std::endl<<std::endl;
        std::cerr << "Command:    quant     Quantitative Trait"                                   << std::endl;
        std::cerr << "            binary    Binary Trait"                                         << std::endl;
        //std::cerr << "            risk-qt   Risk Prediction on continuous traits"                 << std::endl;
        //std::cerr << "            risk-bin  Risk Prediction on dichotomize traits"                << std::endl;
        std::cerr << std::endl;
    }
    else if(mode.compare("version")==0 || mode.compare("--version")==0){
        //Get the version information
        std::cerr << "Version: "<< m_version << std::endl<<std::endl;
    }
    else{
        std::cerr << std::endl;
        std::cerr << "unrecognized command: " << mode << std::endl<< std::endl;
        std::cerr << "Usage: ./SHREK <command> [options]"                                         << std::endl<<std::endl;
        std::cerr << "Command:    quant     Quantitative Trait"                                   << std::endl;
        std::cerr << "            binary    Binary Trait"                                         << std::endl;
        //std::cerr << "            risk-qt   Risk Prediction on continuous traits"                 << std::endl;
       // std::cerr << "            risk-bin  Risk Prediction on dichotomize traits"                << std::endl;
        std::cerr << std::endl;
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
    // Not very efficient here, but should be ok consider the number of possible header
    // Maybe someone who is interested in optimization can optimize this
    bool found = false;
    for(size_t i = 0; i < index.size(); ++i){
        for(size_t j = 0; j < token.size(); ++j){
            if(index[i].compare(token[j])==0){
                found = true;
                indexResult.push_back(j+1); //Don't use it as index, therefore 0 = uninitialized
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
    static const char *optString = "a:R:k:p:b:c:r:l:s:o:f:t:S:N:i:I:X:x:d:L:vunh?";
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
        {"sign",required_argument, NULL,'S'},
        {"nullSign",required_argument, NULL,'N'},
        {"info", required_argument, NULL, 'i'},
        {"ref",required_argument, NULL, 'X' }, //these are for the validation, optional
        {"alt", required_argument, NULL, 'x'}, //these are for the validation, optional
        //General parameters
        {"out",required_argument, NULL,'o'},
        {"impute", required_argument, NULL, 'I'},
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
			case 'c': //chr
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('c');
                break;
			case 'r': //rsid
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('r');
				break;
			case 'l': //loc
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('l');
				break;
			case 'i': //info
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('i');
				break;
            case 'X': //ref
                colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('X');
                break;
            case 'x': //alt
                colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('x');
                break;

            case 's': //statistic
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('s');
                break;
            case 'S': //sign
                colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('S');
                m_signGiven = true;
        //General parameters
			case 'o':
				m_outputPrefix = optarg;
				break;
			case 'f': //maf threshold
                m_maf = atof(optarg);
                m_providedMaf= true;
                break;
			case 'I': //impute threshold
                m_imputeThreshold = atof(optarg);
                break;
			case 'N': //The null value for sign of statistic
                m_nullSign = atof(optarg);
                break;
			case 't': //The number of threads
				m_thread = atoi(optarg);
				break;
			case 'd': //The block distance
				m_distance = atoi(optarg);
				break;
			case 'L': //The list of regions
				m_regionList = optarg;
                break;
			case 'v': //Whether if we should validate the SNPs
				m_validate = true;
				break;
			case 'u': //Whether if the input is p-value or summary statistics
                m_isPvalue = true;
                break;
			case 'n': //Whether to conduct the LD correction
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

    //The problem here is that the sign of the summary statistic can be the same column as the summary statistics
    std::vector<size_t> index;
    getIndex(colNamesForIndex, m_pValueFileName, index);
    if(index.size() != colNamesForIndex.size()){
        throw "Some field(s) not found in header, please check";
    }

    // For all index, most of them should not be a duplicate of each other.
    std::map<size_t, char> dupCheck;
    bool duplicated = false;
    for(size_t i = 0; i < typeForIndex.size(); ++i){
        switch(typeForIndex[i]){
            case 'c':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_chrIndex = index[i];
                    dupCheck[index[i]] = 'c';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'l':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_bpIndex = index[i];
                    dupCheck[index[i]] = 'l';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'r':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_rsIndex = index[i];
                    dupCheck[index[i]] = 'r';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'i':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_imputeInfo = index[i];
                    dupCheck[index[i]] = 'i';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'X':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_ref = index[i];
                    dupCheck[index[i]] = 'X';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'x':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_alt = index[i];
                    dupCheck[index[i]] = 'x';
                }
                else{
                    duplicated = true;
                }
                break;
            case 's':
                if(dupCheck.find(index[i])==dupCheck.end() || dupCheck[index[i]] == 'S'){
                    m_stats = index[i];
                    dupCheck[index[i]] = 's';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'S':
                if(dupCheck.find(index[i])==dupCheck.end() || dupCheck[index[i]] == 's'){
                    m_signIndex = index[i];
                    m_signGiven =true;
                    dupCheck[index[i]] = 'S';
                }
                else{
                    duplicated = true;
                }
                break;
        }
    }

    if(duplicated){
        throw "Duplicated fields, please check your input";
    }

    // Now perform sanity check
    bool error= false;
    std::string errorLog = "";
    if(m_caseSize ==0){
        error = true;
        errorLog.append("Number of case cannot be 0\n");
    }
    if(m_controlSize==0){
        error = true;
        errorLog.append("Number of control cannot be 0\n");
    }
    if(!m_providedPrevalence){
        error = true;
        errorLog.append("The population prevalence is required\n");
    }
    else if(m_prevalence <0 || m_prevalence> 1){
        error =true;
        errorLog.append("The population prevalence must be between 0 or 1\n");
    }
    if(m_pValueFileName.empty()){
        error = true;
        errorLog.append("The p-value file is required\n");
    }
    if(m_ldFilePrefix.empty()){
        error = true;
        errorLog.append("The reference panel is required\n");
    }
    if(m_validate){
        if(m_chrIndex==0){
            error=true;
            errorLog.append("Column for chromosome information not provided\n");
        }
        if(m_rsIndex==0){
            error=true;
            errorLog.append("Column for rs ID not provided\n");
        }
        if(m_bpIndex==0){
            error = true;
            errorLog.append("Column for SNP coordinate not provided\n");
        }
        if(m_stats==0){
            error = true;
            errorLog.append("Column for statistics not provided\n");
        }
        if(m_imputeInfo!= 0 && (m_imputeThreshold < 0 ||m_imputeThreshold> 1) ){
            error = true;
            errorLog.append("Imputation INFO must be between 0 and 1\n");
        }
        if(m_providedMaf && (m_maf < 0 || m_maf > 1)){
            error = true;
            errorLog.append("MAF must be between 0 and 1\n");
        }
    }
    else if(m_stats==0){
        error = true;
        errorLog.append("Summary statistics must be provided\n");
    }
    if(error){
        throw errorLog;
    }

}

void Command::quantitativeProcess(int argc, char* argv[]){
    static const char *optString = "e:a:A:b:p:c:r:l:s:S:N:i:I:X:x:d:f:L:o:t:nuvh?";
	static const struct option longOpts[]={
	    //Qt specific parameter
		{"extreme",required_argument,NULL,'e'},
        {"sampleSize",required_argument,NULL,'a'},
        {"sampleIndex",required_argument,NULL,'A'},
        //File parameters
		{"bfile",required_argument,NULL,'b'},
		{"pfile",required_argument,NULL,'p'},
		{"chr",required_argument,NULL,'c'},
		{"rs",required_argument,NULL,'r'},
		{"loc",required_argument,NULL,'l'},
		{"stats", required_argument, NULL, 's'},
		{"sign",required_argument,NULL,'S'},
        {"nullSign",required_argument, NULL,'N'},
        {"info", required_argument, NULL, 'i'},
        {"impute", required_argument, NULL, 'I'},
        {"ref",required_argument, NULL, 'X' }, //these are for the validation, optional
        {"alt", required_argument, NULL, 'x'}, //these are for the validation, optional
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
			case 'a':
				m_sampleSize= atoi(optarg);
				m_provideSampleSize = true;
				break;
			case 'A':
				//m_sampleSizeIndex = atoi(optarg)-1;
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('A');
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
                break;
			case 'r':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('r');
				break;
			case 'l':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('l');
				break;
            case 'S':
                colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('S');
                m_signGiven = true;
            case 's':
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('s');
                break;
			case 'i': //info
			    colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('i');
				break;
            case 'X': //ref
                colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('X');
                break;
            case 'x': //alt
                colNamesForIndex.push_back(optarg);
                typeForIndex.push_back('x');
                break;
            case 'N':
                m_nullSign = atof(optarg);
                break;
			case 'I': //impute threshold
                m_imputeThreshold = atof(optarg);
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
    std::map<size_t, char> dupCheck;
    bool duplicated = false;
    for(size_t i = 0; i < typeForIndex.size(); ++i){
        switch(typeForIndex[i]){
            case 'A':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_sampleSizeIndex = index[i];
                    dupCheck[index[i]] = 'A';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'c':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_chrIndex = index[i];
                    dupCheck[index[i]] = 'c';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'l':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_bpIndex = index[i];
                    dupCheck[index[i]] = 'l';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'r':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_rsIndex = index[i];
                    dupCheck[index[i]] = 'r';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'i':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_imputeInfo = index[i];
                    dupCheck[index[i]] = 'i';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'X':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_ref = index[i];
                    dupCheck[index[i]] = 'X';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'x':
                if(dupCheck.find(index[i])==dupCheck.end()){
                    m_alt = index[i];
                    dupCheck[index[i]] = 'x';
                }
                else{
                    duplicated = true;
                }
                break;
            case 's':
                if(dupCheck.find(index[i])==dupCheck.end() || dupCheck[index[i]] == 'S'){
                    m_stats = index[i];
                    dupCheck[index[i]] = 's';
                }
                else{
                    duplicated = true;
                }
                break;
            case 'S':
                if(dupCheck.find(index[i])==dupCheck.end() || dupCheck[index[i]] == 's'){
                    m_signIndex = index[i];
                    m_signGiven =true;
                    dupCheck[index[i]] = 'S';
                }
                else{
                    duplicated = true;
                }
                break;
        }
    }

    if(duplicated){
        throw "Duplicated fields, please check your input";
    }

    //Now perform input sanity check;
    bool error= false;
    std::string errorLog = "";

    if(m_provideExtremeAdjustment){
        if(m_extremeAdjust < 0){
            error = true;
            errorLog.append("Extreme adjustment must be bigger than 0\n");
        }
    }
    if(m_provideSampleSize){
        if(m_sampleSize==0){
            error = true;
            errorLog.append("Sample size must be bigger than 0\n");
        }
    }
    else{
        if(m_sampleSizeIndex==0){
            error = true;
            errorLog.append("Sample size information must be provided\n");
        }
    }
    if(m_pValueFileName.empty()){
        error = true;
        errorLog.append("The p-value file is required\n");
    }
    if(m_ldFilePrefix.empty()){
        error = true;
        errorLog.append("The reference panel is required\n");
    }
    if(m_validate){
        if(m_chrIndex==0){
            error=true;
            errorLog.append("Column for chromosome information not provided\n");
        }
        if(m_rsIndex==0){
            error=true;
            errorLog.append("Column for rs ID not provided\n");
        }
        if(m_bpIndex==0){
            error = true;
            errorLog.append("Column for SNP coordinate not provided\n");
        }
        if(m_stats==0){
            error = true;
            errorLog.append("Column for statistics not provided\n");
        }
        if(m_imputeInfo!= 0 && (m_imputeThreshold < 0 ||m_imputeThreshold> 1) ){
            error = true;
            errorLog.append("Imputation INFO must be between 0 and 1\n");
        }
        if(m_providedMaf && (m_maf < 0 || m_maf > 1)){
            error = true;
            errorLog.append("MAF must be between 0 and 1\n");
        }
    }
    else if(m_stats==0){
        error = true;
        errorLog.append("Summary statistics must be provided\n");
    }
    if(error){
        throw errorLog;
    }

}

/*
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

*/


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


