#include "command.h"

Command::Command(){}

Command::~Command(){}
void Command::usage(){
    fprintf(stderr, "Usage:   shrek <command> [options]\n\n");
    fprintf(stderr, "Command: quant        Quantitative Trait\n");
    fprintf(stderr, "         binary       Binary Trait\n");
}

void Command::generalOptions(){

    fprintf(stderr, "Options:\n");
    fprintf(stderr, "         -o | --output      Output prefix [stdout]\n");
    fprintf(stderr, "         -L | --region      Region information, Format: \n");
    fprintf(stderr, "                            name_1:bed_file_1, name2:bed_file_2,...\n");
    //Although we prefer true, we will set default to false to give a clear parameter information
    fprintf(stderr, "         -b | --block       Size of Block [%u]\n", m_sizeOfBlock);
    fprintf(stderr, "         -e | --correct     Perform LD correction [False]\n");
    fprintf(stderr, "         -k | --keep        Keep ambiguous SNPs (e.g. G|C in reference panel and C|G in\n");
    fprintf(stderr, "                            p-value file [False]\n");
    fprintf(stderr, "         -u | --pvalue      P-values instead of summary statistics are provided [False]\n");
    fprintf(stderr, "         -t | --thread      Number of thread used [%u]\n", m_nThread);
    fprintf(stderr, "         -f | --maf         MAF threshold for reference SNPs\n");
    fprintf(stderr, "         -I | --impute      Imputation info score threshold[%f]\n", m_infoThreshold);
    fprintf(stderr, "         -h | --help        Print this help\n");
}

void Command::btUsage(){
    fprintf(stderr, "Usage:   shrek binary [options] \n\n");
    fprintf(stderr, "Required Arguments:\n");
    fprintf(stderr, "         -p | --pfile       P-value file name\n");
    fprintf(stderr, "         -r | --bfile       Reference Panel file prefix\n"); //Currently we only support plink format
    fprintf(stderr, "         -K | --prevalence  Population Prevalence of the trait\n");
    fprintf(stderr, "         -s | --stat        Column name for summary statistic / p-value\n");
    fprintf(stderr, "         -S | --sign        Column name for direction of effect\n");
    fprintf(stderr, "                            e.g. OR, Z e.t.c\n");
    fprintf(stderr, "         -U | --null        Null for direction of effect [0]\n");
    fprintf(stderr, "                            e.g. -U 1 for OR, 0 for others\n");
    fprintf(stderr, "         -v | --conIndex    Column containing the number of controls\n");
    fprintf(stderr, "         -V | --nControl    Number of controls\n");
    fprintf(stderr, "                            Only functional when -v not provided\n");
    fprintf(stderr, "         -w | --caseIndex   Column containing the number of cases\n");
    fprintf(stderr, "         -W | --nCase       Number of cases\n");
    fprintf(stderr, "                            Only functional when -w not provided\n");
    fprintf(stderr, "         -c | --chr         Column name for chromosome \n");
    fprintf(stderr, "         -m | --rs          Column name for rsID \n");
    fprintf(stderr, "         -l | --loc         Column name for coordinate \n");
    fprintf(stderr, "         -a | --ref         Column name for reference allele\n");
    fprintf(stderr, "         -A | --alt         Column name for alternative allele\n");
    fprintf(stderr, "         -i | --info        Column name for impute info score\n");
    generalOptions();
}
void Command::qtUsage(){
    fprintf(stderr, "Usage:   shrek quant  [options] \n\n");
    fprintf(stderr, "Required Arguments:\n");
    fprintf(stderr, "         -p | --pfile       P-value file name\n");
    fprintf(stderr, "         -r | --bFile       Reference Panel file prefix\n"); //Currently we only support plink format
    fprintf(stderr, "         -x | --extreme     Extreme Adjustment, can be calculated as: \n");
    fprintf(stderr, "                            Variance before selection / Variance after selection\n");
    fprintf(stderr, "         -n | --sampIndex   Column name for the number of sample\n");
    fprintf(stderr, "         -N | --nSample     Number of sample in the study\n");
    fprintf(stderr, "                            Only functional when -n not provided\n");
    fprintf(stderr, "         -s | --stat        Column name for summary statistic / p-value\n");
    fprintf(stderr, "         -S | --sign        Column name for direction of effect\n");
    fprintf(stderr, "                            e.g. OR, Z e.t.c\n");
    fprintf(stderr, "         -U | --null        Null for direction of effect [0]\n");
    fprintf(stderr, "                            e.g. -U 1 for OR, 0 for others\n");
    fprintf(stderr, "         -c | --chr         Column name for chromosome \n");
    fprintf(stderr, "         -m | --rs          Column name for rsID \n");
    fprintf(stderr, "         -l | --loc         Column name for coordinate \n");
    fprintf(stderr, "         -A | --ref         Column name for reference allele\n");
    fprintf(stderr, "         -a | --alt         Column name for alternative allele\n");
    fprintf(stderr, "         -i | --info        Column name for impute info score\n");
    generalOptions();
}

//This will return whether we should continue our work
bool Command::parseCommand(int argc, char *argv[]){
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: shrek (Tool for SNP heritability estimation)\n");
    fprintf(stderr, "Version: %f\n\n", m_version);
    if(argc <= 1) usage();
    else{
        return processCode(argc, argv);
    }
    return false;
}


// Return whether if we can continue onwards
// True = Can Progress
// False = Stop programme
bool Command::processCode(int argc,char *argv[]){
    //Do both case control and quantitative parsing together consider how extensive their overlaps are
    std::string mode =argv[1];
    if(argc==2){
        if(mode.compare("quant")==0) qtUsage();
        else if(mode.compare("binary")==0) btUsage();
        else{
            //Incorrect mode
            fprintf(stderr, "Undefined mode %s\n\n", mode.c_str());
            usage();
        }
        return false;
    }
    if(mode.compare("quant")==0) m_qt=true;
    else if(mode.compare("binary")==0) m_caseControl = true;
    static const char *optString = "o:L:b:ekut:f:I:p:r:K:x:s:S:v:V:w:W:c:m:l:a:A:i:n:N:U:h?";
    static const struct option longOpts[]={
	    //Qt specific parameter
		{"output",required_argument,NULL,'o'},
		{"region",required_argument,NULL,'L'},
		{"block",required_argument,NULL,'b'},
		{"correct",no_argument,NULL,'e'},
		{"keep",no_argument,NULL,'k'},
		{"pvalue",no_argument,NULL,'u'},
		{"thread",required_argument,NULL,'t'},
		{"maf",required_argument,NULL,'f'},
		{"impute",required_argument,NULL,'I'},
		{"pfile",required_argument,NULL,'p'},
		{"bfile",required_argument,NULL,'r'},
		{"prevalence",required_argument,NULL,'K'},
		{"extreme",required_argument,NULL,'x'},
		{"stat",required_argument,NULL,'s'},
		{"sign",required_argument,NULL,'S'},
		{"conIndex",required_argument,NULL,'v'},
		{"nControl",required_argument,NULL,'V'},
		{"caseIndex",required_argument,NULL,'w'},
		{"nCase",required_argument,NULL,'W'},
		{"chr",required_argument,NULL,'c'},
		{"rs",required_argument,NULL,'m'},
		{"loc",required_argument,NULL,'l'},
		{"ref",required_argument,NULL,'A'},
		{"alt",required_argument,NULL,'a'},
		{"info",required_argument,NULL,'i'},
		{"sampIndex",required_argument,NULL,'n'},
		{"nSample",required_argument,NULL,'N'},
		{"null",required_argument,NULL,'U'},
		{"help",no_argument,NULL,'h'},
		{NULL, 0, 0, 0}
	};

    bool error = false;
    int longIndex=0;
	int opt = 0;
	opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
	std::vector<char> type;
    std::vector<std::string> columnName;
    std::map<std::string, char> duplication;
	//Start reading all the parameters and perform the qc at the same time
    while(opt!=-1){
		switch(opt){
            case 'o':
                m_outputPrefix = optarg;
                break;
            case 'L':
                m_regionList = optarg;
                break;
            case 'b':
                m_sizeOfBlock = atoi(optarg);
                if(m_sizeOfBlock < 1){
                    error = true;
                    fprintf(stderr, "Size of block must be greater than 0: %s\n", optarg);
                }
                break;
            case 'e':
                m_ldCorrect = true;
                break;
            case 'k':
                m_keepAmbiguous = true;
                break;
            case 'u':
                m_isPvalue = true;
                break;
            case 't':
                m_nThread = atoi(optarg);
                if(m_nThread < 1){
                    error = true;
                    fprintf(stderr, "Number of thread must be greater than 0: %s\n", optarg);
                }
                break;
            case 'f':
                m_mafThreshold = atof(optarg);
                if(m_mafThreshold < 0.0 || m_mafThreshold > 1.0){
                    error = true;
                    fprintf(stderr, "MAF threshold must be between 0 and 1: %s\n", optarg);
                }
                break;
            case 'I':
                m_infoThreshold = atof(optarg);
                if(m_infoThreshold < 0.0 || m_infoThreshold > 1.0){
                    error = true;
                    fprintf(stderr, "INFO threshold must be between 0 and 1: %s\n", optarg);
                }
                break;
            case 'p':
                m_pValueFileName = optarg;
                if(m_pValueFileName.empty() || !usefulTools::fileExists(m_pValueFileName)){
                    error = true;
                    fprintf(stderr, "P-value file is required: %s\n", optarg);
                }
                break;
            case 'r':
                m_referenceFilePrefix = optarg;
                if(m_referenceFilePrefix.empty()){
                    error = true;
                    fprintf(stderr, "Reference panel must be provided: %s\n", optarg);
                }
                else if(!usefulTools::fileExists(m_referenceFilePrefix+".bed")){
                    error = true;
                    fprintf(stderr, "bed file not found: %s.bed\n", optarg);
                }
                else if(!usefulTools::fileExists(m_referenceFilePrefix+".fam")){
                    error = true;
                    fprintf(stderr, "fam file not found: %s.fam\n", optarg);

                }
                else if(!usefulTools::fileExists(m_referenceFilePrefix+".bim")){
                    error = true;
                    fprintf(stderr, "bim file not found: %s.bim\n", optarg);

                }
                break;
            case 'K':
                if(mode.compare("quant")==0){
                    error=true;
                    fprintf(stderr, "Population prevalence not required for quantitative trait analysis\n");
                }
                else if(mode.compare("binary")==0){
                    m_prevalence = atof(optarg);
                    if(m_prevalence < 0.0 || m_prevalence >1.0){
                        error = true;
                        fprintf(stderr, "Population prevalence must be between 0 and 1: %s\n", optarg);
                    }
                }
                break;
            case 'x':
                if(mode.compare("binary")==0){
                    error=true;
                    fprintf(stderr, "Extreme adjustment not applicable for binary trait analysis\n");
                }
                else if(mode.compare("quant")==0){
                    m_extremeAdjust = atof(optarg);
                    if(m_extremeAdjust< 0.0){
                        error = true;
                        fprintf(stderr, "Extreme adjustment value must be larger than 0: %f\n", m_extremeAdjust);
                    }
                }
                break;
            case 'V':
                if(mode.compare("quant")==0){
                    error=true;
                    fprintf(stderr, "Number of controls is only used for binary trait analysis\n");
                    fprintf(stderr, "Please use -N for quantitative trait analysis\n");
                }
                else if(mode.compare("binary")==0){
                    m_nControl = atoi(optarg);
                    if(m_nControl <=0){
                        error = true;
                        fprintf(stderr, "Number of controls must be larger than 0: %i", m_nControl);
                    }
                }
                break;
            case 'W':
                if(mode.compare("quant")==0){
                    error=true;
                    fprintf(stderr, "Number of cases is only used for binary trait analysis\n");
                    fprintf(stderr, "Please use -N for quantitative trait analysis\n");
                }
                else if(mode.compare("binary")==0){
                    m_nCase = atoi(optarg);
                    if(m_nCase <=0){
                        error = true;
                        fprintf(stderr, "Number of cases must be larger than 0: %i\n", m_nCase);
                    }
                }
                break;
            case 'N':
                if(mode.compare("binary")==0){
                    error=true;
                    fprintf(stderr, "Number of samples is only used for quantitative trait analysis\n");
                    fprintf(stderr, "Please use -V and -W for binary trait analysis\n");
                }
                else if(mode.compare("quant")==0){
                    m_nSample = atoi(optarg);
                    if(m_nSample <=0){
                        error = true;
                        fprintf(stderr, "Number of samples must be larger than 0: %i\n", m_nSample);
                    }
                }
                break;
            case 'U':
                m_signNull = atof(optarg);
                break;
            //Problem for the index based input is that they can be provided BEFORE the p-value file is provided
            //Therefore we must do the QC post hoc
            case 's':
                if(duplication.find(optarg)==duplication.end() || duplication[optarg]=='S'){
                    //This is ok, because Sign and statistic can be on the same column
                    duplication[optarg] = 's';
                    type.push_back('s');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'S':
                if(duplication.find(optarg)==duplication.end() || duplication[optarg]=='s'){
                    duplication[optarg] = 'S';
                    type.push_back('S');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'v':
                if(duplication.find(optarg)==duplication.end()){
                    duplication[optarg]='v';
                    type.push_back('v');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'w':
                if(duplication.find(optarg)==duplication.end()){
                    duplication[optarg]='w';
                    type.push_back('w');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'c':
                if(duplication.find(optarg)==duplication.end()){
                    duplication[optarg]='c';
                    type.push_back('c');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'm':
                if(duplication.find(optarg)==duplication.end()){
                    duplication[optarg]='m';
                    type.push_back('m');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'l':
                if(duplication.find(optarg)==duplication.end()){
                    duplication[optarg]='l';
                    type.push_back('l');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'a':
                if(duplication.find(optarg)==duplication.end()){
                    duplication[optarg]='a';
                    type.push_back('a');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'A':
                if(duplication.find(optarg)==duplication.end()){
                    duplication[optarg]='A';
                    type.push_back('A');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'i':
                if(duplication.find(optarg)==duplication.end()){
                    duplication[optarg]='i';
                    type.push_back('i');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
            case 'n':
                if(duplication.find(optarg)==duplication.end()){
                    duplication[optarg]='n';
                    type.push_back('n');
                    columnName.push_back(optarg);
                }
                else{
                    error = true;
                    fprintf(stderr, "Duplicated column name found, please check your input: %s\n", optarg);
                }
                break;
    		case 'h':
			case '?':
			    if(mode.compare("quant")==0) qtUsage();
                else if(mode.compare("binary")==0) btUsage();
                return false; //This is not an error, just tell the programme to end
				break;
			default:
				throw "Undefined operator, please use --help for more information!";
		}
		opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    }

    // now all the basic QC has been completed
    // The only QC left is to make sure all the columns are here
    std::ifstream pvalueFile;
    pvalueFile.open(m_pValueFileName.c_str());
    if(!pvalueFile.is_open()){
        error = false;
        fprintf(stderr, "Cannot open p-value file %s\n", m_pValueFileName.c_str());
    }
    else{
        std::string line;
        std::getline(pvalueFile, line);
        pvalueFile.close();
        std::vector<std::string> token;
        usefulTools::tokenizer(line, "\t ", &token); //Assume the header is separated by tab  or space
        // Therefore the name of the headers cannot contain space

        //Now obtain the index of the corresponding columns
        bool identified = false;
        for(size_t i = 0; i < columnName.size(); ++i){
            for(size_t j = 0; j < token.size(); ++j){
                if(columnName[i].compare(token[j])==0){
                    (j+1 > m_maxIndex)? m_maxIndex=j+1: m_maxIndex=m_maxIndex; //Again, reserve 0 (we are going to use it as bound anyway)
                    identified = true;
                    switch(type[i]){
                        case 's':
                            m_summaryStatisticIndex = j+1;
                            break;
                        case 'S':
                            m_signIndex = j+1;
                            break;
                        case 'v':
                            m_controlIndex=j+1;
                            break;
                        case 'w':
                            m_caseIndex=j+1;
                            break;
                        case 'c':
                            m_chrIndex = j+1;
                            break;
                        case 'm':
                            m_rsIndex = j+1;
                            break;
                        case 'l':
                            m_bpIndex = j+1;
                            break;
                        case 'a':
                            m_altIndex=j+1;
                            break;
                        case 'A':
                            m_refIndex=j+1;
                            break;
                        case 'i':
                            m_imputeInfoIndex=j+1;
                            break;
                        case 'n':
                            m_nSampleIndex=j+1;
                            break;
                        default:
                            assert(false && "Undefined type of parameter");
                    }
                    break;
                }
            }
            if(!identified){
                error = true;
                fprintf(stderr, "Cannot find column: %s\n", columnName[i].c_str());
            }
            identified = false;
        }
        //Add the default
        for(size_t i=0; i < token.size(); ++i){
            if(duplication.find(token[i])==duplication.end()){
                if(token[i].compare("A1")==0){
                    fprintf(stderr, "Reference allele set to default: A1\n");
                    m_refIndex=i+1;
                    m_maxIndex= (i+1 > m_maxIndex)? i+1: m_maxIndex;
                }
                else if (token[i].compare("A2")==0){
                    fprintf(stderr, "Alternative allele set to default: A2\n");
                    m_altIndex=i+1;
                    m_maxIndex= (i+1 > m_maxIndex)? i+1: m_maxIndex;
                }
                else if (token[i].compare("T")==0 && m_summaryStatisticIndex==0 && m_qt){
                    fprintf(stderr, "Test statistic detected: T\n");
                    m_summaryStatisticIndex=i+1;
                    m_maxIndex= (i+1 > m_maxIndex)? i+1: m_maxIndex;
                    if(m_signIndex==0){
                        fprintf(stderr, "Will also use as sign\n");
                        m_signIndex = i+1;
                    }
                }
                else if(token[i].compare("CHR")==0){
                    fprintf(stderr, "Chromosome header set to default: CHR\n");
                    m_chrIndex= i+1;
                    m_maxIndex= (i+1 > m_maxIndex)? i+1: m_maxIndex;
                }
                else if (token[i].compare("BP")==0){
                    fprintf(stderr, "Coordinate header set to default: BP\n");
                    m_bpIndex = i+1;
                    m_maxIndex= (i+1 > m_maxIndex)? i+1: m_maxIndex;
                }
                else if (token[i].compare("SNP")==0){
                    fprintf(stderr, "SNP ID header set to default: SNP\n");
                    m_rsIndex = i+1;
                    m_maxIndex= (i+1 > m_maxIndex)? i+1: m_maxIndex;
                }
                else if(m_qt && m_summaryStatisticIndex==0  && token[i].compare("P")==0 ){
                    fprintf(stderr, "P-value detected but not provided: P\n");
                    fprintf(stderr, "Will use default value\n");
                    m_summaryStatisticIndex= i+1;
                    m_maxIndex= (i+1 > m_maxIndex)? i+1: m_maxIndex;
                    m_isPvalue=true;
                }
            }
        }
    }
    // Check whether if all the required inputs are provided
    if(m_pValueFileName.empty()){
        error = true;
        fprintf(stderr, "P-value file is required\n");
    }
    if(m_referenceFilePrefix.empty()){
        error = true;
        fprintf(stderr, "Reference panel must be provided\n");
    }
    if(m_chrIndex==0){
        error = true;
        fprintf(stderr, "Column of chromosome must be provided\n");
    }
    if(m_rsIndex==0){
        error = true;
        fprintf(stderr, "Column of rsID must be provided\n");
    }
    if(m_bpIndex==0){
        error = true;
        fprintf(stderr, "Column of coordinates must be provided\n");
    }
    if(m_refIndex==0 || m_altIndex==0){
        fprintf(stderr, "WARNING: Column for reference/alternative allele not provided.  SNPs validation can\n");
        fprintf(stderr, "         only be performed based on coordinates which can introduce errors\n");
    }
    if(m_imputeInfoIndex==0){
        fprintf(stderr, "WARNING: Column for impute info not provided. Will not perform filtering\n");
    }
    if(m_summaryStatisticIndex==0){
        error = true;
        fprintf(stderr, "Column name for summary statistic / p-value must be provided\n");
    }
    if(m_signIndex==0){
        fprintf(stderr, "WARNING: Column for sign not provided. Variance estimation might not be accurate\n");
    }
    if(m_caseControl){
        if(m_nCase < 1 && m_caseIndex==0){
            error = true;
            fprintf(stderr, "Number of cases must be provided\n");
        }
        else if(m_caseIndex != 0 && m_nCase >= 1){
                m_nCase = 0;
            fprintf(stderr, "Input number of cases not used as column for case number is provided\n");
        }
        if(m_nControl < 1 && m_controlIndex==0){
            error = true;
            fprintf(stderr, "Number of controls must be provided\n");
        }
        else if(m_controlIndex != 0 && m_nControl >= 1){
            m_nControl = 0;
            fprintf(stderr, "Input number of controls not used as column for control number is provided\n");
        }
        else if(m_prevalence < 0){
            error = true;
            fprintf(stderr, "Population prevalence must be provided\n");
        }
    }
    else if(m_qt){
        if(m_nSampleIndex == 0 && m_nSample < 1){
            error = true;
            fprintf(stderr, "Sample size must be provided\n");
        }
        else if(m_nSampleIndex != 0 && m_nSample >= 1){
            m_nSample = 0;
            fprintf(stderr, "Input sample size not used as column for sample size is provided\n");
        }

    }

    //Check if the imputation filtering is ok
    if(m_imputeInfoIndex!=0){
        if(m_infoThreshold==1.0){
            fprintf(stderr, "All SNPs will be filtered as maximum INFO score can only be 1\n");
            error=true;
        }
        else if(m_infoThreshold==0.0){
            fprintf(stderr, "INFO threshold not provided (-I), will set to default: 0.8\n");
            m_imputeInfoIndex=0.8;
        }
    }
    // If there is no error then everything is ok
    if(error){
        throw "Problem with parameter input(s)";
    }
    return !error;
}
