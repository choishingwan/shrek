#include "genotype.h"

size_t Genotype::m_sampleNum;

Genotype::Genotype(){
    m_bitSize = sizeof(mlong)*CHAR_BIT;
    m_requiredBit = Genotype::m_sampleNum*2;
    m_genotype = new mlong [(m_requiredBit /(m_bitSize))+1];
    m_missing = new mlong [(m_requiredBit /(m_bitSize))+1];
    memset(m_genotype, 0x0,((m_requiredBit /(m_bitSize))+1)*sizeof(mlong));
    memset(m_missing, 0x0,((m_requiredBit /(m_bitSize))+1)*sizeof(mlong));
}


Genotype::~Genotype(){
    delete [] m_genotype;
    delete [] m_missing;
}

void Genotype::GetbothR(const Genotype &snpB, const bool correction, double &r, double &rSq) const {
    size_t range = (m_requiredBit /(m_bitSize))+1;
    mlong acc1 = 0;
    mlong acc2 = 0;
    mlong acc11 = 0;
    mlong acc22 = 0;
    mlong acc = 0;
    double res[5];
    res[0] = (double)m_sampleNum;
    res[1] = -(double)snpB.m_nonMissSample;
    res[2] = -(double)m_nonMissSample;
    res[3] = res[1];
    res[4] = res[2];
    double n = 0.0;
    for(size_t i = 0; i < range; ++i){
        mlong loader1 = m_genotype[i];
        mlong loader2 = snpB.m_genotype[i];
        mlong sum1 = snpB.m_missing[i];
        mlong sum2 = m_missing[i];
        n+= (double)__builtin_popcountll(sum1&sum2)/2.0;
        mlong sum12 = (loader1 | loader2) & m1;
        sum1 = sum1&loader1;
        sum2 = sum2&loader2;
        mlong sum11 = sum1 & m1;
        mlong sum22 = sum2 & m1;
        loader1= ((~(m1+sum12)) & (loader1^loader2));
        sum12 = sum12|loader1;
        sum1 = ((sum1 & m2) + ((sum1 >> 2) & m2));
        sum2 = ((sum2 & m2) + ((sum2 >> 2) & m2));
        sum12 = ((sum12 & m2) + ((sum12 >> 2) & m2));
        sum11 =((sum11 & m2) + ((sum11 >> 2) & m2));
        sum22 =((sum22 & m2) + ((sum22 >> 2) & m2));
        acc1 += ((sum1&m4)+((sum1>>4)& m4));
        acc2 += ((sum2& m4) + ((sum2 >> 4) & m4));
        acc11 += ((sum11& m4)+ ((sum11 >> 4) & m4));
        acc22 += ((sum22& m4)+ ((sum22 >> 4) & m4));
        acc += ((sum12& m4)+ ((sum12 >> 4) & m4));
    }
    acc1 = (acc1&m8)+((acc1>>8)&m8);
    acc2 = (acc2&m8)+((acc2>>8)&m8);
    acc = (acc&m8)+((acc>>8)&m8);
    acc11 = ((acc11+(acc11>>8)) & m8);
    acc22 = ((acc22+(acc22>>8)) & m8);
    res[0] -= (acc*0x1000100010001LLU >> 48);
    res[1] += (acc1*0x1000100010001LLU >> 48);
    res[2] += (acc2*0x1000100010001LLU >> 48);
    res[3] += (acc11*0x1000100010001LLU >> 48);
    res[4] += (acc22*0x1000100010001LLU >> 48);
    double dxx = res[1];
    double dyy = res[2];
    //double n = (double)m_nonMissSample - ((double)m_sampleNum-(double)snpB.m_nonMissSample);
    double cov12 = res[0] * n - dxx * dyy;
    dxx = (res[3] * n + dxx * dxx) * (res[4] * n + dyy * dyy);
    if(dxx==0.0){
        r=0.0;
        rSq=0.0;
    }
    else{
        r=cov12 / sqrt(dxx);
        rSq =(cov12 * cov12) / dxx;
        std::cerr << "Checking: " << rSq << std::endl;
        if(correction){
            r= r*(1+(1-r*r)/(2*(n-4))); //POPA
            rSq=rSq-1.0/(2.0*n); //Weir & Hill
        }
    }
}


void Genotype::SetsampleNum(size_t sampleNum){
    if(sampleNum < 2){
        throw "ERROR! Sample number must be bigger than 1!";
    }
    Genotype::m_sampleNum = sampleNum;
}



void Genotype::AddsampleGenotype(const int first, const int second, const size_t sampleIndex){
    if(first==second && first==1){ //hom ref 10
        m_missing[(sampleIndex*2)/(m_bitSize)] = m_missing[(sampleIndex*2)/(m_bitSize)]  | 0x3LLU << ((sampleIndex*2)% (m_bitSize));
        m_genotype[(sampleIndex*2)/(m_bitSize)] = m_genotype[(sampleIndex*2)/(m_bitSize)]  | 0x2LLU << ((sampleIndex*2)% (m_bitSize));
        m_nonMissSample++;
    }
    else if(first==second && first==0){//hom alt 00
        m_missing[(sampleIndex*2)/(m_bitSize)] = m_missing[(sampleIndex*2)/(m_bitSize)]  | 0x3LLU << ((sampleIndex*2)% (m_bitSize));
        m_nonMissSample++;
    }
    else if(first==0 && second==1){ //Het 01
        m_missing[(sampleIndex*2)/(m_bitSize)] = m_missing[(sampleIndex*2)/(m_bitSize)]  | 0x3LLU << ((sampleIndex*2)% (m_bitSize));
        m_genotype[(sampleIndex*2)/(m_bitSize)] = m_genotype[(sampleIndex*2)/(m_bitSize)]  | 0x1LLU << ((sampleIndex*2)% (m_bitSize));
        m_nonMissSample++;
    }
    else if(first==1 && second ==0){ //missing
        m_genotype[(sampleIndex*2)/(m_bitSize)] = m_genotype[(sampleIndex*2)/(m_bitSize)]  | 0x1LLU << ((sampleIndex*2)% (m_bitSize));
    }
}

