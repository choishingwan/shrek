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
    mlong loader1, loader2, sum1, sum2, sum11, sum12, sum22;
	uint32_t final_sum1 = 0;
	uint32_t final_sum2 = 0;
	uint32_t final_sum11 = 0;
	uint32_t final_sum22 = 0;
	uint32_t final_sum12 = 0;
	double return_vals[5];
    return_vals[0] = (double)m_sampleNum;
    return_vals[1] = -(double)snpB.m_nonMissSample;
    return_vals[2] = -(double)m_nonMissSample;
    return_vals[3] = return_vals[1];
    return_vals[4] = return_vals[2];
    unsigned int N =0;
    for(size_t i = 0; i < range;){
        loader1 = m_genotype[i];
        loader2 = snpB.m_genotype[i];
        sum1 = snpB.m_missing[i];
        sum2 = m_missing[i];
		i++;
		N+= __builtin_popcountll(sum1&sum2)/2;
		sum12 = (loader1 | loader2) & FIVEMASK;
		sum1 = sum1 & loader1;
		sum2 = sum2 & loader2;
		loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
		sum12 = sum12 | loader1;
		sum11 = sum1 & FIVEMASK;
		sum22 = sum2 & FIVEMASK;
		sum1 = (sum1 & 0x33333333) + ((sum1 >> 2) & 0x33333333);
		sum2 = (sum2 & 0x33333333) + ((sum2 >> 2) & 0x33333333);
		sum12 = (sum12 & 0x33333333) + ((sum12 >> 2) & 0x33333333);
		mlong tmp_sum1=0 , tmp_sum2=0;
		if(i < range){
	    	loader1 = m_genotype[i];
	    	loader2 = snpB.m_genotype[i];
	    	tmp_sum1 = snpB.m_missing[i];
	    	tmp_sum2 =m_missing[i];
			N+= __builtin_popcountll(tmp_sum1&tmp_sum2)/2;
		}
		else{
			loader1 = 0;
			loader2 = 0;
		}
	    i++;
	    mlong tmp_sum12 = (loader1 | loader2) & FIVEMASK;
		tmp_sum1 = tmp_sum1 & loader1;
		tmp_sum2 = tmp_sum2 & loader2;
		loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
		tmp_sum12 = tmp_sum12 | loader1;
		sum11 += tmp_sum1 & FIVEMASK;
		sum22 += tmp_sum2 & FIVEMASK;
		sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
		sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
		sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);
	    if(i < range){
			loader1 = m_genotype[i];
			loader2 = snpB.m_genotype[i];
			tmp_sum1 = snpB.m_missing[i];
			tmp_sum2 =m_missing[i];
			N+= __builtin_popcountll(tmp_sum1&tmp_sum2)/2;
		}
		else{
			loader1=0;
			loader2=0;
			tmp_sum1=0;
			tmp_sum2=0;
		}
		i++;
		tmp_sum12 = (loader1 | loader2) & FIVEMASK;
		tmp_sum1 = tmp_sum1 & loader1;
		tmp_sum2 = tmp_sum2 & loader2;
		loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
		tmp_sum12 = tmp_sum12 | loader1;
		sum11 += tmp_sum1 & FIVEMASK;
		sum22 += tmp_sum2 & FIVEMASK;
		sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
		sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
		sum11 = (sum11 & 0x33333333) + ((sum11 >> 2) & 0x33333333);
		sum22 = (sum22 & 0x33333333) + ((sum22 >> 2) & 0x33333333);
		sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);
		sum1 = (sum1 & 0x0f0f0f0f) + ((sum1 >> 4) & 0x0f0f0f0f);
		sum2 = (sum2 & 0x0f0f0f0f) + ((sum2 >> 4) & 0x0f0f0f0f);
		sum11 = (sum11 & 0x0f0f0f0f) + ((sum11 >> 4) & 0x0f0f0f0f);
		sum22 = (sum22 & 0x0f0f0f0f) + ((sum22 >> 4) & 0x0f0f0f0f);
		sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);
		final_sum1 += (sum1 * 0x01010101) >> 24;
		final_sum2 += (sum2 * 0x01010101) >> 24;
		final_sum11 += (sum11 * 0x01010101) >> 24;
		final_sum22 += (sum22 * 0x01010101) >> 24;
		final_sum12 += (sum12 * 0x01010101) >> 24;
	}

	return_vals[0] -= final_sum12;
	return_vals[1] += final_sum1;
	return_vals[2] += final_sum2;
	return_vals[3] += final_sum11;
	return_vals[4] += final_sum22;

    double dxx = return_vals[1];
    double dyy = return_vals[2];
    double n = N;
    double cov12 = return_vals[0] * n - dxx * dyy;
    dxx = (return_vals[3] * n + dxx * dxx) * (return_vals[4] * n + dyy * dyy);
    if(dxx !=0.0){
        r=cov12 / sqrt(dxx);
        rSq =(cov12 * cov12) / dxx;
        if(correction){
            r= r*(1.0+(1.0-r*r)/(2.0*(n-4.0))); //POPA
            //rSq = 1.0-((m_sampleNum-3.0)/(m_sampleNum-2.0))*(1.0-rSq)*(1.0+(2.0*(1.0-rSq))/(m_sampleNum-3.3)); //OP5 from Shieh
            rSq=rSq-1.0/(2.0*n); //Weir & Hill
    //        rSq=rSq-(1.0-rSq)/(n-2.0); //Weir & Hill
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

