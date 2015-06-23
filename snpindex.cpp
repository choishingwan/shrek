#include "snpindex.h"

SnpIndex::SnpIndex()
{
	m_index = std::map<std::string, size_t>();
	m_isInitialized=false;
}

SnpIndex::~SnpIndex()
{
	//dtor
}

bool SnpIndex::valid() {
    if(!m_isInitialized) throw "SnpIndex isn't initialized";
    return m_iter != m_index.end();
}
bool SnpIndex::contains(std::string input) { return m_index.find(input) != m_index.end(); }
size_t SnpIndex::value() const {
    if(!m_isInitialized) throw "SnpIndex isn't initialized";
    return m_iter->second;
}
size_t SnpIndex::value (std::string key){
	return m_index[key];
}
bool SnpIndex::next(){
    if(!m_isInitialized) throw "SnpIndex isn't initialized";
	m_iter++;
	return m_iter != m_index.end();
}

size_t SnpIndex::size() const { return m_index.size(); }
void SnpIndex::increment(std::string key){ m_index[key]++; }
void SnpIndex::set(std::string key, size_t value){ m_index[key] = value; };
void SnpIndex::init(){
    m_isInitialized=true;
    m_iter = m_index.begin();
}
void SnpIndex::print(){
    std::map<std::string, size_t>::iterator iter;
    for(iter = m_index.begin(); iter != m_index.end(); ++iter){
        std::cerr << iter->first << "\t" << iter->second << std::endl;
    }
}



