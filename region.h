// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef REGION_H
#define REGION_H

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <complex>
#include "interval.h"
#include "usefulTools.h"





/** \class Region
 *  \brief The region class, store the variance information of each region
 *
 *  The main goal of this class is to allow user to provide multiple bed
 *  files when they would like to run the programme so that heritability
 *  estimation can be performed on all the regions at once. e.g. different
 *  snp categories, different gene sets etc.
 *  This class will also form the basis to construct the flag information
 *  of the Snps.
 *  The first region is the default region, which contains all the snps that *  have the LD information.
 */
class Region
{
	public:
	    /** Default constructor */
	    Region();
	    /** Default destructor */
		virtual ~Region();
        /** Based on the given command line input, read the region information
         *  The region should be in the format:
         * \<Name\>:\<Bed File\>,<Name\>:\<Bed File\>...
         */
        void generateRegion(std::string regionList);
        /** Remove all the interval pointers as they are no longer required */
        void clean();
        /** Add variance to the i th region */
        void Addvariance(double const var, size_t i);
        /** Add the additional variance to the i th region */
        void AddadditionVariance(double const addVar, size_t i);
        /** Return the number of regions */
        size_t GetnumRegion() const;
        /** Return the chromosome information of the jth interval in the ith region */
        std::string Getchr(size_t i, size_t j) const;
        /** Return the start coordinate of the jth interval in the ith region */
        size_t Getstart(size_t i, size_t j) const;
        /** Return the last coordinate of the jth interval in the ith region */
        size_t Getend(size_t i, size_t j) const;
        /** Return the size of the ith region */
        size_t GetintervalSize(size_t i) const;
        /** Return the name of the i th region */
        std::string Getname(size_t i) const;
        /** Return the variance of the i th region given the heritability */
        double Getvariance(double heritability, size_t i) const;
	protected:
	private:
        std::vector<std::vector<Interval*> > m_intervalList;
        std::vector<std::string> m_names;
        std::vector<double> m_variance;
        std::vector<double> m_additionVariance;
};

#endif // REGION_H
