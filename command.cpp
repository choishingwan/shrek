#include "command.h"

size_t Command::Getthread() const { return m_thread; }
size_t Command::GetblockSize() const { return m_blockSize; }
size_t Command::Getdistance() const { return m_distance; }
size_t Command::GetstepSize() const { return m_stepSize; }
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
double Command::Getinflation() const { return m_inflation; }
double Command::Getmaf() const { return m_maf; }
bool Command::ldCorrect() const { return m_ldCorrection; }
bool Command::validate() const { return m_validate; }
bool Command::isPvalue() const { return m_isPvalue; }
bool Command::provideSampleSize() const { return m_provideSampleSize; }
bool Command::quantitative() const { return m_quantitative; }
bool Command::caseControl() const { return m_caseControl; }
std::string Command::GetoutputPrefix() const { return m_outputPrefix; }
std::string Command::GetpValueFileName() const { return m_pValueFileName; }
std::string Command::GetgenotypeFilePrefix() const { return m_genotypeFilePrefix; }
std::string Command::GetregionList() const { return m_regionList; }
std::string Command::GetprogrammeName() const { return m_programmeName; }



Command::Command(int argc, char* argv[])
{

    if(argc == 1){
		std::cerr << "You have not provided any arguments. Please provide all the required arguments" << std::endl;
		printBriefUsage();
		exit(-1);
	}
    int threadDefault = 1;
    int blockSizeDefault = 420;
    bool providedPrevalence = false;
    m_ldCorrection = true;
    m_validate = false;
	m_isPvalue = false;
	m_provideSampleSize = false;
    m_quantitative = false;
    m_caseControl = false;
    m_provideSampleSize =false;
    m_thread = threadDefault;
    m_blockSize = blockSizeDefault;
	m_chrIndex = 0;
    m_rsIndex = 1;
    m_bpIndex = 2;
    m_sampleSizeIndex = 3;
	m_maf = 0.0;

    m_genotypeFilePrefix = "";
    m_pValueFileName = "";
	m_regionList="";
	m_programmeName =argv[0];
	static const char *optString = "t:b:d:e:s:a:q:c:r:x:k:m:nvuo:p:g:L:h?";
	static const struct option longOpts[]={
		{"thread", required_argument, NULL, 't'},
        {"blockSize", required_argument, NULL, 'b'},
        {"distance", required_argument, NULL, 'd'},
        {"step", required_argument, NULL, 'e'},
		{"sampleSize", required_argument, NULL, 's'},
		{"case", required_argument, NULL, 0},
		{"control", required_argument, NULL, 0},
        {"cc", required_argument, NULL, 'a'},
        {"quant", required_argument, NULL, 'q'},
		{"bp", required_argument, NULL, 0},
        {"chr", required_argument, NULL, 'c'},
		{"rsid", required_argument, NULL, 'r'},
		{"sampleIndex", required_argument, NULL, 'x'},
		{"prevalence", required_argument, NULL, 'k'},
		{"maf", required_argument, NULL, 'm'},
		{"no_correct", no_argument, NULL, 'n'},
		{"validate", no_argument, NULL, 'v'},
		{"pvalue", no_argument, NULL, 'u'},
		{"out", required_argument, NULL, 'o'},
		{"pfile", required_argument, NULL, 'p'},
		{"genotype", required_argument, NULL, 'g'},
		{"region", required_argument, NULL, 'L'},

        {"help", no_argument, NULL, 'h'},
		{NULL, 0, 0, 0}
	};

	bool error = false;
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
				m_blockSize = atoi(optarg);
				break;
			case 'd':
				m_distance = atoi(optarg);
				break;
			case 'e':
				m_stepSize = atoi(optarg);
				break;
            case 's':
				m_sampleSize= atoi(optarg);
				break;
            case 'a':
				m_cIndex= atoi(optarg)-1;
				m_caseControl = true;
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
				m_provideSampleSize = true;
				break;
			case 'k':
				m_prevalence = atof(optarg);
				providedPrevalence =true;
				break;
            case 'm':
                m_maf = atof(optarg);
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
            case 'g':
                m_genotypeFilePrefix = optarg;
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
	if(m_blockSize <= 2){
		std::cerr << "Undefined block size: " << m_blockSize << ". Block size must be greater than 2, set block size to default: " << blockSizeDefault << std::endl;
		m_blockSize = blockSizeDefault;
	}
	else if(m_blockSize %3 != 0){
        std::cerr << "Block size should be dividable by 3 to avoid complex algorithm. Block size now changed to " ;
        m_blockSize = m_blockSize-(m_blockSize%3);
        std::cerr << m_blockSize << std::endl;
	}
    if(m_caseControl && m_quantitative){
		std::cerr << "You must select either case control study or quantitative traits, but not both or none" << std::endl;
		error=true;
    }
	else if(m_rsIndex <= 0 || m_chrIndex <0 || m_bpIndex <0 || m_sampleSizeIndex <0){
        std::cerr << "The column number must be a positive integer" << std::endl;
        error=true;
	}
    else if(m_caseControl && (m_cIndex==m_chrIndex || m_cIndex == m_rsIndex || m_cIndex==m_bpIndex)){
        std::cerr << "The column should be unique!" << std::endl;
        error=true;
    }
    else if(m_quantitative && (m_tIndex==m_chrIndex || m_tIndex == m_rsIndex || m_tIndex==m_bpIndex || (m_provideSampleSize && m_tIndex==m_sampleSizeIndex))){
        std::cerr << "The column should be unique!" << std::endl;
        error=true;
    }
    if(m_caseControl && (m_caseSize == 0 || m_controlSize == 0 )){
		std::cerr << "Require to provide number of case and control for case control study" << std::endl;
		std::cerr << "And the number of case and number of control must both be larger than 0" << std::endl;
		error=true;
	}
	else if(m_caseControl && !providedPrevalence){
		std::cerr << "Require to provide the prevalence for case control study" << std::endl;
		error = true;
	}
	else if(m_caseControl){
		if(m_caseSize + m_controlSize == 0){
			std::cerr << "The Sum of Case and control are equal to zero! Please check that your input is correct" << std::endl;
			error=true;
		}
	}
	if(m_chrIndex == m_rsIndex || m_chrIndex == m_bpIndex || m_rsIndex == m_bpIndex){
        std::cerr << "The number of columns must be unique!" << std::endl;
        error=true;
	}
	else if(m_provideSampleSize && (m_sampleSizeIndex == m_rsIndex || m_sampleSizeIndex == m_bpIndex || m_sampleSizeIndex == m_chrIndex) ){
        std::cerr << "The number of columns must be unique!" << std::endl;
        error=true;
	}
	if(m_genotypeFilePrefix.empty()){
		std::cerr << "You must provide the genotype file for LD construction!" << std::endl;
		error = true;
	}
	else if(!usefulTools::fileExists(m_genotypeFilePrefix+".bed") &&
            !usefulTools::fileExists(m_genotypeFilePrefix+".bim") &&
            !usefulTools::fileExists(m_genotypeFilePrefix+".fam")){
        std::cerr << "Cannot open the genotype file(s), please check if they exists" << std::endl;
        std::cerr << "Genotype prefix: " << m_genotypeFilePrefix << std::endl;
        error = true;
	}
	if(m_pValueFileName.empty()){
        error=true;
        std::cerr << "You must provide the p-value file for the analysis" << std::endl;
    }
    else if(!usefulTools::fileExists(m_pValueFileName)){
        error = true;
        std::cerr << "Cannot open p-value file: " << m_pValueFileName << std::endl;
        std::cerr << "Please check if they exists" << std::endl;

    }
    if(m_maf < 0.0 || m_maf > 1.0){
        error = true;
        std::cerr << "Invalid maf! maf must be within the range of 0.0 to 1.0. Your input: " << m_maf << std::endl;
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

}

void Command::printUsage(){

    std::cerr << "------------------------------------------------------------------------------"  << std::endl;
    std::cerr << "| Heritability Estimate using summary statistic                               |" << std::endl;
    std::cerr << "| version 0.01                                                                |" << std::endl;
    std::cerr << "| (C) 2014 Johnny Kwan, Sam Choi                                              |" << std::endl;
    std::cerr << "| The University of Hong Kong                                                 |" << std::endl;
    std::cerr << "| Haven't figure out which license                                            |" << std::endl;
    std::cerr << "------------------------------------------------------------------------------"  << std::endl;
    std::cerr << "usage: ./Jest [-p <p-value_file>] [--tstat <t-stat_index> | --chi <chi_index>]"  << std::endl;
    std::cerr << "              [ --g <genotype_prefix> | --l <linkage_file> ] ..."                << std::endl;
    std::cerr << "Required options: "                                                              << std::endl;
    std::cerr << "  -p,--pfile       The p-value file. The first 3 fields of the file must be"     << std::endl;
    std::cerr << "                   <Chr>    <rsID>    <BP>"                                      << std::endl;
    std::cerr << "                   and test statistic for each snp must also be provided"        << std::endl;
    std::cerr << "  -c,--chr         The column number of chromosome in the p-value file"          << std::endl;
    std::cerr << "  --rsid           The column number of rsid in the p-value file"                << std::endl;
    std::cerr << "  --bp             The column number of rsid coordinate in the p-value file"     << std::endl;
    std::cerr << "  -g,--genotype    The genotype file prefix. The programme will use this to "    << std::endl;
    std::cerr << "                   calculate the LD matrix. Will require the fam, bim and bed"   << std::endl;
    std::cerr << "                   file. Please try to perform quality control beforehand "      << std::endl;
    std::cerr << "  -n,--no_correct  Turn of LD correction. The LD correction is used when "       << std::endl;
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
    std::cerr << "  -b,--blockSize   The size of the block of analysis. As it is computationally " << std::endl;
    std::cerr << "                   difficult to perform decomposition for the whole genome, "    << std::endl;
    std::cerr << "                   we break each chromosomes into a small blocks and perform a"  << std::endl;
    std::cerr << "                   sliding window approach. Please note that this option will "  << std::endl;
    std::cerr << "                   be the main determinant of memory usage and run time and "    << std::endl;
    std::cerr << "                   will have a effect to the result. If the block size is too "  << std::endl;
    std::cerr << "                   small, the result will tend to be over-estimated; if the "    << std::endl;
    std::cerr << "                   block size is too big, it might take forever to finish the "  << std::endl;
    std::cerr << "                   analysis. As we use a step size of blockSize/3, we require "  << std::endl;
    std::cerr << "                   the blockSize to be divisible by 3 to avoid off-by-one "      << std::endl;
    std::cerr << "                   error. Shall the used provide value not divisible by 3, we "  << std::endl;
    std::cerr << "                   will change it to a value divisible by 3."                    << std::endl;
    std::cerr << "  -d,--distance    The flanking distance of each snps to be included. Run time " << std::endl;
    std::cerr << "                   increase exponentially with the selection of this number. "   << std::endl;
    std::cerr << "                   This number will also affect the prediction of the programme" << std::endl;
    std::cerr << "  -e,--step        The step size of the window. The larger it is, the faster  "  << std::endl;
    std::cerr << "                   the programme runs. "                                         << std::endl;
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
    std::cerr 	<< "Jest\tCoding by Sam CHOI\tMethod by Johnny KWAN" <<std::endl
        << "===============================================================" << std::endl
		<< "Performing analysis using the following parameters: " << std::endl
		<< "===============================================================" << std::endl
		<< "Essential Input  " <<std::endl;
		std::cerr	<< "Genotype File Prefix : " << m_genotypeFilePrefix << std::endl;
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
				<< "Number of Thread     : " << m_thread << std::endl
				<< "Block Size           : " << m_blockSize << std::endl;
    if(m_ldCorrection){
        std::cerr << "Use LD correction    : True" << std::endl;
    }
    else{
        std::cerr << "Use LD correction    : False" << std::endl;
    }
	std::cerr	<< "Number of regions    : " << regionMessage << std::endl;
	std::cerr << std::endl << std::endl;

}

