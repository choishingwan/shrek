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
#include <boost/ptr_container/ptr_vector.hpp>
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
        Region();
        virtual ~Region();
        void generateRegion(std::string regionList);
        void clean();
        size_t getNumRegion() const {return m_names.size();};
        inline std::string getChr(size_t i, size_t j) const {return m_intervalList.at(i).at(j).Getchr();};
        inline size_t getStart(size_t i, size_t j) const {return m_intervalList.at(i).at(j).Getstart();};
        inline size_t getEnd(size_t i, size_t j) const {return m_intervalList.at(i).at(j).Getend();};
        inline size_t getIntervalSize(size_t i) const {return m_intervalList.at(i).size();};
        inline std::string getName(size_t i) const {return m_names.at(i);};
        protected:
	private:
        std::vector<boost::ptr_vector<Interval> > m_intervalList;
        std::vector<std::string> m_names;
};

#endif // REGION_H
