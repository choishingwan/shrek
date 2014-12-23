#include "snpestimation.h"

SnpEstimation::SnpEstimation(GenotypeFileHandler *genotypeFileHandler, SnpIndex *snpIndex, std::vector<Snp*> *snpList):m_genotypeFileHandler(genotypeFileHandler), m_snpIndex(snpIndex), m_snpList(snpList){};

SnpEstimation::~SnpEstimation()
{
	//dtor
}
