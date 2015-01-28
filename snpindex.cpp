#include "snpindex.h"

SnpIndex::SnpIndex()
{
	m_index = std::map<std::string, size_t>();
}

SnpIndex::~SnpIndex()
{
	//dtor
}

bool SnpIndex::valid() { return m_iter != m_index.end(); }
bool SnpIndex::find(std::string input) { return m_index.find(input) != m_index.end(); }
size_t SnpIndex::value() const { return m_iter->second; }
size_t SnpIndex::value (std::string key){
	return m_index[key];
}
bool SnpIndex::next(){
	m_iter++;
	return m_iter != m_index.end();
}

size_t SnpIndex::size() const { return m_index.size(); }
void SnpIndex::increment(std::string key){ m_index[key]++; }
void SnpIndex::set(std::string key, size_t value){ m_index[key] = value; };
void SnpIndex::init(){ m_iter = m_index.begin(); }
void SnpIndex::print(){
    std::map<std::string, size_t>::iterator iter;
    for(iter = m_index.begin(); iter != m_index.end(); ++iter){
        std::cerr << iter->first << "\t" << iter->second << std::endl;
    }
}



