#include "decomposition.h"

Decomposition::Decomposition(SnpIndex *snpIndex, std::vector<Snp*> *snpList, Linkage *linkageMatrix, size_t thread):m_snpIndex(snpIndex), m_snpList(snpList), m_linkage(linkageMatrix), m_thread(thread){}

Decomposition::~Decomposition()
{
	//dtor
}
