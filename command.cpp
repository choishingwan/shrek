#include "command.h"

size_t Command::Getthread() const { return m_thread; }
size_t Command::GetminBlock() const { return m_minBlock; }
size_t Command::GetmaxBlock() const { return m_maxBlock; }
size_t Command::GetsampleSize() const { return m_sampleSize; }
size_t Command::GetcaseSize() const { return m_caseSize; }
size_t Command::GetcontrolSize() const { return m_controlSize; }
size_t Command::GetcIndex() const { return m_cIndex; }
size_t Command::GettIndex() const { return m_tIndex; }
size_t Command::GetbpIndex() const { return m_bpIndex; }
size_t Command::GetchrIndex() const { return m_chrIndex; }
size_t Command::GetrsIndex() const { return m_rsIndex; }
size_t Command::GetsampleSizeIndex() const { return m_sampleSizeIndex; }
double Command::Getprevalence() const { return m_prevalence; }
double Command::Getmaf() const { return m_maf; }
double Command::GetextremeAdjust() const { return m_extremeAdjust; }
bool Command::ldCorrect() const { return m_ldCorrection; }
bool Command::validate() const { return m_validate; }
bool Command::isPvalue() const { return m_isPvalue; }
bool Command::provideSampleSize() const { return m_provideSampleSize; }
bool Command::quantitative() const { return m_quantitative; }
bool Command::caseControl() const { return m_caseControl; }
bool Command::maxBlockSet() const { return m_maxBlockSet; }
std::string Command::GetoutputPrefix() const { return m_outputPrefix; }
std::string Command::GetpValueFileName() const { return m_pValueFileName; }
std::string Command::GetldFilePrefix() const { return m_ldFilePrefix; }
std::string Command::GetregionList() const { return m_regionList; }
std::string Command::GetprogrammeName() const { return m_programmeName; }
std::string Command::GetdirectionFile() const { return m_directionFile; }



Command::Command(int argc, char* argv[], bool &error)
{
    error = false;
    if(argc == 1){
		std::cerr << "You have not provided any arguments. Please provide all the required arguments" << std::endl;
		printBriefUsage();
		error = true;
		return;
	}
    int threadDefault = 1;
    bool providedPrevalence = false;
    bool providedMaf= false;
    bool provideExtremeAdjustment = false;
    m_ldCorrection = true;
    m_validate = false;
	m_isPvalue = false;
	m_provideSampleSize = false;
    m_quantitative = false;
    m_caseControl = false;
    m_maxBlockSet = false;
    m_thread = threadDefault;
    m_chrIndex = 0;
    m_rsIndex = 1;
    m_bpIndex = 2;
    m_sampleSizeIndex = 3;
	m_maf = -1.0;
    m_maxBlock = 0;
    m_minBlock = 0;
    m_ldFilePrefix = "";
    m_pValueFileName = "";
	m_regionList="";
	m_programmeName =argv[0];
    m_sampleSize=0;
	m_caseSize=0;
	m_controlSize=0;
	m_cIndex=7;
    m_tIndex=7;
    m_prevalence=1.0;
    m_extremeAdjust=1.0;
    m_outputPrefix="";
	m_pValueFileName="";
    m_ldFilePrefix="";
	m_regionList="";
	m_directionFile="";

	static const char *optString = "t:b:B:s:a:q:c:r:x:k:e:m:d:nvuo:p:l:L:h?";
	static const struct option longOpts[]={
		{"thread", required_argument, NULL, 't'},
        {"minBlock", required_argument, NULL, 'b'},
        {"maxBlock", required_argument, NULL, 'B'},
		{"sampleSize", required_argument, NULL, 's'},
		{"case", required_argument, NULL, 0},
		{"control", required_argument, NULL, 0},
        {"cc", required_argument, NULL, 'a'},
        {"quant", required_argument, NULL, 'q'},
		{"bp", required_argument, NULL, 0},
        {"chr", required_argument, NULL, 'c'},
		{"rs", required_argument, NULL, 'r'},
		{"sampleIndex", required_argument, NULL, 'x'},
		{"prevalence", required_argument, NULL, 'k'},
		{"maf", required_argument, NULL, 'm'},
		{"no_correct", no_argument, NULL, 'n'},
		{"validate", no_argument, NULL, 'v'},
		{"pvalue", no_argument, NULL, 'u'},
		{"direction", no_argument, NULL, 'u'},
		{"out", required_argument, NULL, 'o'},
		{"pfile", required_argument, NULL, 'p'},
		{"linkage", required_argument, NULL, 'l'},
		{"region", required_argument, NULL, 'L'},
		{"extreme", required_argument, NULL, 'e'},
        {"help", no_argument, NULL, 'h'},
		{NULL, 0, 0, 0}
	};

    int longIndex;
	int opt = 0;
	std::string interpret="";
	opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    while(opt!=-1){
		switch(opt){
			case 0:
				interpret = longOpts[longIndex].name;
				if(interpret.compare("control") == 0){
					m_controlSize= atoi(optarg);
				}
				if(interpret.compare("case") == 0){
					m_caseSize= atoi(optarg);
				}
				if(interpret.compare("bp")==0){
                    m_bpIndex=atoi(optarg)-1;
				}
				break;
			case 't':
				m_thread = atoi(optarg);
				break;
			case 'b':
				m_minBlock = atoi(optarg);
				break;
			case 'B':
				m_maxBlock = atoi(optarg);
				m_maxBlockSet =true;
				break;
			case 's':
				m_sampleSize= atoi(optarg);
				m_provideSampleSize = true;
				break;
            case 'a':
				m_cIndex= atoi(optarg)-1;
				m_caseControl = true;
				break;
            case 'd':
				m_directionFile = optarg;
				break;
            case 'e':
				m_extremeAdjust = atof(optarg);
				provideExtremeAdjustment = true;
				break;
            case 'q':
				m_tIndex= atoi(optarg)-1;
				m_quantitative = true;
				break;
            case 'c':
                m_chrIndex = atoi(optarg)-1;
                break;
			case 'r':
				m_rsIndex = atoi(optarg)-1;
				break;
			case 'x':
				m_sampleSizeIndex = atoi(optarg)-1;
				break;
			case 'k':
				m_prevalence = atof(optarg);
				providedPrevalence =true;
				break;
            case 'm':
                m_maf = atof(optarg);
                providedMaf= true;
                break;
			case 'n':
				m_ldCorrection = false;
				break;
			case 'v':
				m_validate = true;
				break;
            case 'u':
                m_isPvalue = true;
                break;
			case 'o':
				m_outputPrefix = optarg;
				break;
			case 'p':
				m_pValueFileName = optarg;
				break;
            case 'l':
                m_ldFilePrefix = optarg;
                break;
            case 'L':
				m_regionList = optarg;
                break;
			case 'h':
			case '?':
 				printUsage();
				exit(0);
				break;
			default:
				std::cerr << "Undefined operator, please use --help for more information!" << std::endl;
				std::cerr << optarg << std::endl;
				printBriefUsage();
				exit(-1);
		}
		opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    }

	if(m_thread <= 0){
		std::cerr << "Undefined number of thread(s): " << m_thread << ". Number of thread(s) must be greater than 0, set number of thread(s) to default: " << threadDefault << std::endl;
		m_thread = threadDefault;
	}

    if(m_maxBlockSet && m_minBlock > m_maxBlock){ //Enabled max block, so we need to have minBlock smaller than maxBlock
        error = true;
        std::cerr << "Maximum block size enabled. Therefore minimum block size must be smaller than or equal to maximum block size" << std::endl;
        std::cerr << "Minimum block size: " << m_minBlock << std::endl;
        std::cerr << "Maximum block size: " << m_maxBlock << std::endl;
    }

    if(m_maxBlockSet && m_maxBlock == 0){
		error = true;
        std::cerr << "Maximum block size enabled. The maximum block size must be bigger than 0" << std::endl;
    }
    if(m_caseControl && m_caseSize == 0){
        error = true;
        std::cerr << "Your case control study has 0 case and we cannot perform the analysis on such type of study." << std::endl;
        std::cerr << "Please check your input is correct. Sorry." << std::endl;
    }
    else if(m_caseControl && m_controlSize == 0){
        std::cerr << "Warning! Your case control study has 0 control and we are uncertain how will this affect the result." << std::endl;
        std::cerr << "Please be cautious with the result" << std::endl;
    }
    if(m_quantitative && m_caseControl){
        error = true;
        std::cerr << "You may specify the study as either quantitative or case control study but not both." << std::endl;
        std::cerr << "Please make sure you have the correct input" << std::endl;
    }
    if(m_quantitative && m_provideSampleSize && m_sampleSize <= 0){
        error = true;
        std::cerr << "Sample size provided for the quantitative study is less than or equal to zero" << std::endl;
        std::cerr << "Please check you have the correct input" << std::endl;
    }
    if(m_pValueFileName.empty()){
        error = true;
        std::cerr << "You must provide the p-value file input!" << std::endl;
    }
    else if(!usefulTools::fileExists(m_pValueFileName)){
        error = true;
        std::cerr << "Cannot open the p-value file, please check that the file exists" << std::endl;
    }

    if(m_quantitative && (m_tIndex == m_bpIndex || m_tIndex == m_chrIndex || m_tIndex == m_rsIndex || m_tIndex == m_sampleSizeIndex ||
       m_bpIndex == m_chrIndex || m_bpIndex == m_rsIndex || m_bpIndex == m_sampleSizeIndex ||
       m_chrIndex == m_rsIndex || m_chrIndex == m_sampleSizeIndex ||
       m_rsIndex == m_sampleSizeIndex)){
        error = true;
        std::cerr << "Duplicated index! Please make sure the index are not duplicated!" << std::endl;
        std::cerr << "Statistic index: " << m_tIndex << std::endl;
        std::cerr << "bp index: " << m_bpIndex << std::endl;
        std::cerr << "chr index: " << m_chrIndex << std::endl;
        std::cerr << "rsId index: " << m_rsIndex << std::endl;
        std::cerr << "sample size index: " << m_sampleSizeIndex << std::endl;

    }
    else if(m_caseControl && (m_cIndex == m_bpIndex || m_cIndex == m_chrIndex || m_cIndex == m_rsIndex ||
                              m_bpIndex==m_chrIndex || m_bpIndex == m_rsIndex ||
                              m_chrIndex == m_rsIndex)){
        error = true;
        std::cerr << "Duplicated index! Please make sure the index are not duplicated!" << std::endl;
        std::cerr << "Statistic index: " << m_cIndex << std::endl;
        std::cerr << "bp index: " << m_bpIndex << std::endl;
        std::cerr << "chr index: " << m_chrIndex << std::endl;
        std::cerr << "rsId index: " << m_rsIndex << std::endl;
    }
    if(m_caseControl && !providedPrevalence ){
        error = true;
        std::cerr << "You must provide the prevalence for case control study." << std::endl;
    }
    if(providedMaf && (m_maf < 0.0 || m_maf > 1.0)){
        error = true;
        std::cerr << "maf must be between 0.0 and 1.0" << std::endl;
        std::cerr << "maf input: " << m_maf << std::endl;
    }

    if(m_ldFilePrefix.empty()){
        error = true;
        std::cerr << "Genotype files must be provided for ld calculation" << std::endl;
    }
    if(!error && m_maxBlockSet && m_maxBlock%3 != 0){
        std::cerr << "We prefer a blockSize that can be divided by 3. Will change the max block size to " << m_maxBlock+3-m_maxBlock%3 << std::endl;
        m_maxBlock = m_maxBlock+3-m_maxBlock%3;
    }
	if(!error && m_minBlock > 0 && m_minBlock%3 != 0){
        std::cerr << "We prefer a blockSize that can be divided by 3. Will change the min block size to " << m_minBlock+3-m_minBlock%3 << std::endl;
        m_minBlock = m_minBlock+3-m_minBlock%3;
    }
	if(m_extremeAdjust <= 0.0){
        error = true;
        std::cerr << "The extreme adjustment value should always be bigger than 0" << std::endl;
	}
	else if(provideExtremeAdjustment && m_caseControl){
        std::cerr << "Currently there is no support for extreme phenotype in case control study. Extreme adjustment value will have no effect" << std::endl;
	}
	if(error){
		std::cerr << "Type " << argv[0] << " -h for more information" << std::endl;
		exit(-1);
	}
}

Command::~Command()
{
    //dtor
}

void Command::printBriefUsage(){
    std::cerr << "------------------------------------------------------------------------------"  << std::endl;
    std::cerr << "| Snp HeRitability Estimation Kit                                             |" << std::endl;
    std::cerr << "| version 0.01                                                                |" << std::endl;
    std::cerr << "| (C) 2014 Johnny Kwan, Sam Choi                                              |" << std::endl;
    std::cerr << "| The University of Hong Kong                                                 |" << std::endl;
    std::cerr << "| Haven't figure out which license                                            |" << std::endl;
    std::cerr << "------------------------------------------------------------------------------"  << std::endl;
    std::cerr << "usage: ./Jest [-p <p-value_file>] [--tstat <t-stat_index> | --chi <chi_index>]"  << std::endl;
    std::cerr << "              [ --g <linkage_file_prefix> ] ..."                << std::endl;
    std::cerr << "Required options: "                                                              << std::endl;
    std::cerr << "  -p,--pfile       The p-value file.                              [ Required ]"  << std::endl;
    std::cerr << "  -c,--chr         The column number of chromosome              [ Default: 1 ]"  << std::endl;
    std::cerr << "  --rs             The column number of rsid                    [ Default: 2 ]"  << std::endl;
    std::cerr << "  --bp             The column number of coordinate              [ Default: 3 ]"  << std::endl;
    std::cerr << "  -l,--linkage     The linkage  file prefix.                      [ Required ]"  << std::endl;
    std::cerr << "  -n,--no_correct  Turn off LD correction. "                                     << std::endl;
    std::cerr << "  -m,--maf         The minor allele frequency filtering.      [ Default: off ]"  << std::endl;
    std::cerr << "                                                                              "  << std::endl;
    std::cerr << "Quantitative trait analysis: "                                                   << std::endl;
    std::cerr << "  --quant          The column of statistic. Quantitative trait    [ Required ]"  << std::endl;
    std::cerr << "  -s,--sampleSize  The number of sample used."                                   << std::endl;
    std::cerr << "  -e,--extreme     The extreme phenotype adjustment ratio."                      << std::endl;
    std::cerr << "  --sampleIndex    The sample column index.                     [ Default: 4 ]"  << std::endl;
    std::cerr << "                                                                              "  << std::endl;
    std::cerr << "Case Control analysis: "                                                         << std::endl;
    std::cerr << "  --cc             The column of statistic. Case Control study    [ Required ]"  << std::endl;
    std::cerr << "  -k,--prevalence  Prevalence of the phenotype                    [ Required ]"  << std::endl;
    std::cerr << "  --case           The number of case used in the study           [ Required ]"  << std::endl;
    std::cerr << "  --control        The number of control used in the study        [ Required ]"  << std::endl;
    std::cerr << "                                                                              "  << std::endl;
    std::cerr << "General options: " << std::endl;
    std::cerr << "  -u,--pvalue      Input is p-value "                                            << std::endl;
    std::cerr << "  -d,--direction   File containing the direction of effect "                     << std::endl;
    std::cerr << "  -b,--minBlock    The minimum block size                    [ Default: None ]"  << std::endl;
    std::cerr << "  -B,--maxBlock    The maximum block size                    [ Default: None ]"  << std::endl;
    std::cerr << "  -L,--region      Region information"                                           << std::endl;
    std::cerr << "  -t,--thread      The number of thread                         [ Default: 1 ]"  << std::endl;
    std::cerr << "  -o,--out         The output file prefix."                                      << std::endl;
    std::cerr << "  -h,-?,--help     Display the detail help message"                              << std::endl;
}

void Command::printUsage(){

    std::cerr << "------------------------------------------------------------------------------"  << std::endl;
    std::cerr << "| Snp HeRitability Estimation Kit                                             |" << std::endl;
    std::cerr << "| version 0.01                                                                |" << std::endl;
    std::cerr << "| (C) 2014 Johnny Kwan, Sam Choi                                              |" << std::endl;
    std::cerr << "| The University of Hong Kong                                                 |" << std::endl;
    std::cerr << "| Haven't figure out which license                                            |" << std::endl;
    std::cerr << "------------------------------------------------------------------------------"  << std::endl;
    std::cerr << "usage: ./Jest [-p <p-value_file>] [--tstat <t-stat_index> | --chi <chi_index>]"  << std::endl;
    std::cerr << "              [ --g <genotype_prefix> ] ..."                << std::endl;
    std::cerr << "Required options: "                                                              << std::endl;
    std::cerr << "  -p,--pfile       The p-value file. The first 3 fields of the file must be"     << std::endl;
    std::cerr << "                   <Chr>    <rsID>    <BP>"                                      << std::endl;
    std::cerr << "                   and test statistic for each snp must also be provided"        << std::endl;
    std::cerr << "  -c,--chr         The column number of chromosome in the p-value file"          << std::endl;
    std::cerr << "  --rs             The column number of rsid in the p-value file"                << std::endl;
    std::cerr << "  --bp             The column number of rsid coordinate in the p-value file"     << std::endl;
    std::cerr << "  -l,--linkage     The linkage  file prefix. The programme will use this to "    << std::endl;
    std::cerr << "                   calculate the LD matrix. Will require the fam, bim and bed"   << std::endl;
    std::cerr << "                   file. Please try to perform quality control beforehand "      << std::endl;
    std::cerr << "  -n,--no_correct  Turn off LD correction. The LD correction is used when "      << std::endl;
    std::cerr << "                   using the genotype file to calculate the LD matrix. The "     << std::endl;
    std::cerr << "                   LD correction is performed to adjust for number of sample"    << std::endl;
    std::cerr << "                   used for calculating the LD. It is not recommended to turn"   << std::endl;
    std::cerr << "                   off this function"                                            << std::endl;
    std::cerr << "  -m,--maf         The minor allele frequency filtering. If you provide the "    << std::endl;
    std::cerr << "                   genotype file, you may use only snps passing maf threshold"   << std::endl;
    std::cerr << "                   for the LD calculation. Has no effect otherwise."             << std::endl;
    std::cerr << "                                                                               " << std::endl;
    std::cerr << "Quantitative trait analysis: "                                                   << std::endl;
    std::cerr << "  --quant          The input is quantitative trait. The number should indicate"  << std::endl;
    std::cerr << "                   the field where the statistic or p-value located within the " << std::endl;
    std::cerr << "                   p-value file. "                                               << std::endl;
    std::cerr << "                   For example, for p-value file with the following format: "    << std::endl;
    std::cerr << "                   <Chr>    <rsID>    <BP>    <Statistic> "                      << std::endl;
    std::cerr << "                   you should use --quant 4 "                                    << std::endl;
    std::cerr << "  -s,--sampleSize  The number of sample used for the analysis. If not provided " << std::endl;
    std::cerr << "                   the programme will try to obtain it directly from the "       << std::endl;
    std::cerr << "                   p-value file based on the sample column index information"    << std::endl;
    std::cerr << "  -e,--extreme     The extreme phenotype adjustment value. Should be: "          << std::endl;
    std::cerr << "                   Variance after selection / Variance before selection"         << std::endl;
    std::cerr << "  --sampleIndex    The sample column index. Indicating which column contains "   << std::endl;
    std::cerr << "                   the sample size information "                                 << std::endl;
    std::cerr << "                                                                               " << std::endl;
    std::cerr << "Case Control analysis: "                                                         << std::endl;
    std::cerr << "  --cc             The input is a case control analysis. The user are required"  << std::endl;
    std::cerr << "                   to provide the prevalence of the phenotype and also the "     << std::endl;
    std::cerr << "                   number of case and control used in the study. Similar to "    << std::endl;
    std::cerr << "                   quant, cc indicate location of the statistic or p-value"      << std::endl;
    std::cerr << "                   within the p-value file"                                      << std::endl;
    std::cerr << "  -k,--prevalence  Specify the prevalence of the phenotype"                      << std::endl;
    std::cerr << "  --case           The number of case used in the study"                         << std::endl;
    std::cerr << "  --control        The number of control used in the study"                      << std::endl;
    std::cerr << "                                                                               " << std::endl;
    std::cerr << "General options: " << std::endl;
    std::cerr << "  -u,--pvalue      Indicate whether if the input is p-value or test-statistic  " << std::endl;
    std::cerr << "                   Cannot handle p-value of 0 or 1"                              << std::endl;
    std::cerr << "  -d,--direction   The file containing the direction of significance. If it is " << std::endl;
    std::cerr << "                   not provided, there will be an upward bias of variance. "     << std::endl;
    std::cerr << "                   Format of the file should be: "                               << std::endl;
    std::cerr << "                   <rsID>    <Direction> "                                       << std::endl;
    std::cerr << "                   where the direction should be represented by +1/-1"           << std::endl;
    std::cerr << "  -b,--minBlock    The minimum block size for the sliding window. A large"       << std::endl;
    std::cerr << "                   number will increase the run time exponentially. Window size" << std::endl;
    std::cerr << "                   is the most important factor affecting the run time and "     << std::endl;
    std::cerr << "                   memory usage of the programme. Must be bigger than the "      << std::endl;
    std::cerr << "                   maximum block size if maximum block size is used. Default=0"  << std::endl;
    std::cerr << "  -B,--maxBlock    The maximum block size allowed. This helps to set a limit on" << std::endl;
    std::cerr << "                   the block size used. Will help to avoid using a huge block."  << std::endl;
    std::cerr << "                   However, if the maximum block doesn't cover a full LD block " << std::endl;
    std::cerr << "                   it is likely that the final estimate will be inflated"        << std::endl;
    std::cerr << "                   Default = no maximum block size limit"                        << std::endl;
    std::cerr << "  -L,--region      The region field. You may provide a bed file in the format:"  << std::endl;
    std::cerr << "                   <region name>:<fileName>,<region name>:<fileName>,..."        << std::endl;
    std::cerr << "                   The summary output will provide the per region estimate "     << std::endl;
    std::cerr << "                   where the region are label by their order of input"           << std::endl;
    std::cerr << "  -t,--thread      The number of thread use. The memory requirement of the "     << std::endl;
    std::cerr << "                   programme = thread use x blockSize/3 + (blockSize/3)*2"       << std::endl;
    std::cerr << "  -o,--out         The output file prefix. If provided, the programme will "     << std::endl;
    std::cerr << "                   generate a <output>.sum and <output>.res file where the "     << std::endl;
    std::cerr << "                   <output>.sum file will provide the run summary and the "      << std::endl;
    std::cerr << "                   <output>.res file will provide the per Snp estimate"          << std::endl;
    std::cerr << "  -h,-?,--help     Display the detail help message (This message)"               << std::endl;

}



void Command::printRunSummary(std::string regionMessage){
    std::cerr 	<< "SHREK\tCoding by Sam CHOI\tMethod by Johnny KWAN" <<std::endl
        << "===============================================================" << std::endl
		<< "Performing analysis using the following parameters: " << std::endl
		<< "===============================================================" << std::endl
		<< "Essential Input  " <<std::endl;
		std::cerr	<< "Genotype File Prefix : " << m_ldFilePrefix << std::endl;
	    std::cerr 	<< "P-Value File         : " << m_pValueFileName << std::endl;
    if(m_isPvalue){
	    std::cerr 	<< "Input is P-value     : True" << std::endl;
    }
    else{
        std::cerr 	<< "Input is P-value     : False" << std::endl;
    }
    if(!m_outputPrefix.empty()){
        std::cerr   << "Output File          : " << m_outputPrefix << std::endl;
    }
    if(m_caseControl){
		std::cerr	<< "Mode                 : Case Control " << std::endl
					<< "Number of Case       : " << m_caseSize << std::endl
					<< "Number of Control    : " << m_controlSize << std::endl
					<< "Prevalence           : " << m_prevalence << std::endl << std::endl;
	}
	else{
		std::cerr	<< "Mode                 : Quantitative Trait" << std::endl;
		if(m_provideSampleSize) std::cerr	<< "Sample Size          : " << m_sampleSize << std::endl <<std::endl;
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
	std::cerr	<< "Number of regions    : " << regionMessage << std::endl;
	std::cerr << std::endl << std::endl;

}

