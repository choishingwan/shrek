#include "interval.h"

Interval::Interval(std::string chr, size_t start, size_t end):m_chr(chr), m_start(start), m_end(end){}

Interval::~Interval()
{
    //dtor
}


std::string Interval::Getchr() const { return m_chr; }
size_t Interval::Getstart() const { return m_start; }
size_t Interval::Getend() const { return m_end; }
