#ifndef REGION_H
#define REGION_H

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include "usefulTools.h"



/**
	TODO:
		Re-work on this class and the related stuff. Too many static stuff and I don't like it
**/
class Region
{
	public:
		Region(std::string chr, size_t start, size_t end);
		virtual ~Region();
        static void generateRegion(std::vector<std::vector<Region*> > &regionOut, std::string regionList);
        static void cleanRegion(std::vector<std::vector<Region*> > &regionList);
        static std::vector<std::string> regionNames;
        static std::vector<double> regionVariance;
        static size_t numRegion();
        std::string Getchr() const;
        size_t Getstart() const;
        size_t Getend() const;
	protected:
	private:
        std::string m_chr;
        size_t m_start;
        size_t m_end;

};

#endif // REGION_H
