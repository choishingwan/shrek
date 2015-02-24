#include "genotype.h"

size_t Genotype::m_sampleNum;

Genotype::Genotype(){
	m_bitSize = sizeof(unsigned long long);
	m_requiredBit = Genotype::m_sampleNum*4;
	m_genotypeA = new unsigned long long [(m_requiredBit /(8*m_bitSize))+1];
	m_genotypeB = new unsigned long long [(m_requiredBit /(8*m_bitSize))+1];
	memset(m_genotypeA, 0x0,((m_requiredBit /(8*m_bitSize))+1)*sizeof(unsigned long long));
	memset(m_genotypeB, 0x0,((m_requiredBit /(8*m_bitSize))+1)*sizeof(unsigned long long));
}

Genotype::~Genotype(){
    delete [] m_genotypeA;
    delete [] m_genotypeB;
}

double Genotype::Getr(Genotype* snpB, bool correction){
	size_t range = (m_requiredBit /(8*m_bitSize))+1;
    double rSquare = 0.0;
    size_t numSampleInBlock = 2*m_bitSize;
    size_t i = 0;
	for(; i < range-1;){
		rSquare += (__builtin_popcountll(m_genotypeA[i] & snpB->m_genotypeB[i] )- numSampleInBlock*m_mean*snpB->m_mean)/(m_standardDeviation *snpB->m_standardDeviation);
		i++;
	}
	size_t remainSample = (Genotype::m_sampleNum)%(2*m_bitSize);
	if(remainSample > 0){
		rSquare += (__builtin_popcountll(m_genotypeA[i] & snpB->m_genotypeB[i] )- remainSample*m_mean*snpB->m_mean)/(m_standardDeviation *snpB->m_standardDeviation);
    }
    rSquare *= 1.0/(m_sampleNum-1.0);
	rSquare *= rSquare;
	if(correction){
        return 1.0-((m_sampleNum-3.0)/(m_sampleNum-2.0))*(1.0-rSquare)*(1.0+(2.0*(1.0-rSquare))/(m_sampleNum-3.3));
	}

    return rSquare;
}

void Genotype::SetsampleNum(size_t sampleNum){
	if(sampleNum < 2){
        std::cerr << "ERROR! Sample number must be bigger than 1!" << std::endl;
        exit(-1);
	}
    Genotype::m_sampleNum = sampleNum;
}


void Genotype::Setmean(double mean){ m_mean = mean; }
void Genotype::SetstandardDeviation(double standardDeviation){ m_standardDeviation = standardDeviation; }


void Genotype::AddsampleGenotype(int genotype, size_t sampleIndex){
    switch(genotype){
	case 0:
		break;
	case 1:
		m_genotypeA[(sampleIndex*4)/(8*m_bitSize)] = m_genotypeA[(sampleIndex*4)/(8*m_bitSize)]  | 0x5ull  << ((sampleIndex*4)% (8*m_bitSize));
		m_genotypeB[(sampleIndex*4)/(8*m_bitSize)] = m_genotypeB[(sampleIndex*4)/(8*m_bitSize)]  | 0x3ull << ((sampleIndex*4)% (8*m_bitSize));
		break;
	case 2:
		m_genotypeA[(sampleIndex*4)/(8*m_bitSize)]= m_genotypeA[(sampleIndex*4)/(8*m_bitSize)]  | 0xfull << ((sampleIndex*4)% (8*m_bitSize));
		m_genotypeB[(sampleIndex*4)/(8*m_bitSize)]= m_genotypeB[(sampleIndex*4)/(8*m_bitSize)]  | 0xfull << ((sampleIndex*4)% (8*m_bitSize));
        break;
	default:
        std::cerr << "Undefined genotype: " << genotype << " please check your input!" << std::endl;
        std::cerr << "We expect the bed file to contain only 0, 1 or 2 and nothing else (missing will be considered as 0)" << std::endl;
        exit(-1);
    }
}

void Genotype::clean(std::deque<Genotype*> &genotype, size_t remaining){
	size_t sizeOfGenotype = genotype.size();
	size_t removeCount = sizeOfGenotype-remaining;
    for(size_t i = 0; i < removeCount; ++i){
        Genotype* temp = genotype.front();
        genotype.pop_front();
        delete temp;
    }
}
