// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef SNP_H
#define SNP_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <boost/ptr_container/ptr_vector.hpp>
#include <map>
#include "usefulTools.h"
#include "region.h"
#include "command.h"

/**
 * \class Snp
 * \brief Store Snp information and perform some basic logistic of these information
 */
class Snp
{
public:
        Snp(std::string chr, std::string rs, size_t bp, size_t sampleSize, double original, std::string refAllele, std::string altAllele, int direction);
        virtual ~Snp();
        static void generateSnpList(boost::ptr_vector<Snp> &snpList, const Command &commander);
        static void generateSnpIndex(std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, const Command &commander, const Region &regionList);


        inline std::string getChr() const{return m_chr;};
        inline std::string getRs() const{return m_rs; };
        inline std::string getRef() const{return m_ref; };
        inline std::string getAlt() const{return m_alt; };
        inline size_t getBp() const{return m_bp;};
        inline size_t getSampleSize() const{return m_sampleSize;};
        inline static size_t getMaxSample() {return Snp::MAX_SAMPLE_SIZE; };
        inline int getDirection() const{return m_direction;};
        inline double getBeta() const{ return m_beta; };
        inline double getOriginal() const{ return m_original; };
        inline double getHeritability() const{ return m_heritability; };
        inline double getEffectiveNumber() const{ return m_effectiveNumber; };
        inline double getLDScore() const{ return m_ldscore; };
        inline bool Concordant(std::string chr, size_t bp, std::string rsId) const{  return chr.compare(m_chr) ==0 && bp==m_bp && rsId.compare(m_rs) == 0; }
        inline bool getFlag(size_t i ) const {return m_regionFlag.at(i); };
        inline void setHeritability(double h){ m_heritability = h; };
        inline void setEffectiveNumber(double e){ m_effectiveNumber = e; };
        inline void setLDScore(double s){ m_ldscore = s; };
        Snp(const Snp& that) = delete;
        void setFlag(const size_t i, bool flag);

protected:
private:
        std::string m_chr="";
        std::string m_rs="";
        std::string m_ref="";
        std::string m_alt="";
        size_t m_bp=0;
        size_t m_sampleSize=0;
        static size_t MAX_SAMPLE_SIZE;
        double m_original=0.0;
        int m_direction=1;
        //std::vector<double> m_original;
        //std::vector<double> m_heritability;
        double m_heritability=0.0;
        double m_effectiveNumber=0.0;
        double m_ldscore=0.0;
        //std::vector<double> m_beta;
        double m_beta=0.0;
        //std::vector<bool> m_remove;
        std::vector<bool> m_regionFlag;

        //Functions
        static bool sortSnp(const Snp& i, const Snp& j);
        void computeVarianceExplained(const Command &commander);

};

#endif // SNP_H
