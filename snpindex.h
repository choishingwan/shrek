#ifndef SNPINDEX_H
#define SNPINDEX_H

#include <map>
#include <string>
#include <iostream>

class SnpIndex
{
	public:
		SnpIndex();
		virtual ~SnpIndex();
		size_t value() const;
		size_t value(std::string key);
		size_t size() const;
		void set(std::string key, size_t value);
		void init();
		void increment(std::string key);
		bool find(std::string key);
		bool valid();
		bool next();
		void print();
	protected:
	private:
		std::map<std::string, size_t> m_index;
		std::map<std::string, size_t>::iterator m_iter;
};

#endif // SNPINDEX_H
