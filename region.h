#ifndef REGION_H
#define REGION_H

#include <string>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include "interval.h"
#include "usefulTools.h"

// Each region should only contain the information of their own intervals
// and the corresponding variance
class Region
{
    public:
        /** Default constructor */
        Region(const std::string name);
        Region(const std::string name, const std::string bedFileName);
        /** Default destructor */
        size_t getIntervalSize() const{return m_intervalList.size(); };
        std::string getChr(size_t i) const {return m_intervalList.at(i).getChr();};
        size_t getStart(size_t i) const {return m_intervalList.at(i).getStart();};
        size_t getEnd(size_t i) const {return m_intervalList.at(i).getEnd();};

        virtual ~Region();
        // The aim of this function is to generate the required region list, each region item
        // should contain its corresponding intervals
        static void generateRegionList(boost::ptr_vector<Region> &regionList, const std::string regionInput);
    protected:
    private:
        boost::ptr_vector<Interval> m_intervalList;
        std::string m_name;
        double m_variance=0.0;
        double m_varianceBuff=0.0; //This contain the variance that are not included into the final variance
};

#endif // REGION_H
