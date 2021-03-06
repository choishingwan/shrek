#ifndef SNP_H
#define SNP_H

#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include <fstream>
#include <stdio.h>
#include <errno.h>
#include <stdexcept>
#include <math.h>
#include <algorithm>
#include <map>
#include "usefulTools.h"
#include "region.h"
#include "command.h"

class Snp
{
    public:
        /** Default constructor */
        Snp(std::string chr, std::string rsId, size_t loc, size_t nSample, size_t nCase, size_t nControl, std::string refAllele, std::string altAllele, bool binary, bool quantitative, double statistic, double info, int sign);
        /** Default destructor */
        virtual ~Snp();

        // Static functions
        static void generateSnpList(boost::ptr_vector<Snp> &snpList, const Command &commander);
        static void generateSnpIndex(std::map<std::string, size_t> &snpIndex, boost::ptr_vector<Snp> &snpList, boost::ptr_vector<Region> const &regionList);

        //Getters
        std::string getChr() const {return m_chr; };
        std::string getRs() const{ return m_rsId; };
        size_t getLoc() const{return m_loc;};
        int getNSample() const{return m_nSample;};
        int getNCase() const{return m_nCase;};
        int getNControl() const{return m_nControl;};
        int getSampleSize() const{return m_nCase+m_nControl+m_nSample;};
        std::string getRef() const{return m_ref;};
        std::string getAlt() const{return m_alt;};
        double getTStat() const{return (*m_tStat)/(double)m_tStat.use_count();};
        double getFStat() const{return (*m_fStat)/(double)m_fStat.use_count();};
        double getHeritability() const{ return (*m_heritability)/(double)m_heritability.use_count();};
//        double getHeritability() const{ return m_heritDebug;};
        double getLDSC() const {return (*m_ldScore);};
        double getEffective() const{return m_effectiveNumber;};
        double getInfo() const{return m_infoScore; };
        int getSign() const{return m_sign;};
        bool flag(size_t loc) const {return m_regionFlag.at(loc); };
        char getStatus() const{return m_status;};
        //Setters
        void setFlag(const size_t i, bool flag);
        void setStatus(char status);
        void setHeritability(double herit){ (*m_heritability)=herit;m_heritDebug=herit;};
        void setEffective(double effective){ m_effectiveNumber=effective; };
        //Checking
        // This programme will return whether if the result is concordant and will flip accordingly
        // It will flip the SNP even if it is ambiguous, but will let the caller know through ambig
        bool concordant(const std::string chr, const size_t loc, const std::string rs, std::string refAllele, std::string altAllele, bool &ambig);

        //Perfect LD stuff
        void shareHeritability( Snp& i );

        //Uncertain
        Snp(const Snp& that) = delete; //I honestly can't remember what is this for, I guess this is to disable copying
    protected:
    private:
        //Properties of SNP that we want to know
        std::string m_chr="";
        std::string m_rsId ="";
        size_t m_loc=0;
        int m_nSample=0;
        int m_nCase=0;
        int m_nControl=0;
        std::string m_ref="";
        std::string m_alt="";
        bool m_binary=false;
        bool m_quantitative=false;
        double m_effectiveNumber =0.0; // This is for the calculation of heuristic variance
        double m_infoScore = 0.0;
        int m_sign  = 0; //When sign = 0, it means there is no sign given
        // The status flag, indicates whether if it is:
        // include (I)
        // not found in the reference (R)
        // invalid when compared to reference (V)
        // In perfect LD (L)
        // remove due to MAF (F)
        // Need to write these in the manual
        char m_status='R';
        /** Special members use for managing the perfect LD stuff **/
        std::shared_ptr<double> m_tStat; //Average of the Z score of all SNPs with perfect LD
        std::shared_ptr<double> m_fStat; //Average of the f-state of all SNPs with perfect LD
        std::shared_ptr<double> m_heritability; //The master heritability
        double m_heritDebug = 0.0;
        std::shared_ptr<double> m_ldScore;
        Snp* m_targetClass=this;
        std::vector<bool> m_regionFlag;

        static bool sortSnp(const Snp& i, const Snp& j);
        void computeSummaryStatistic(bool isP);
        void flip(); // Flip the sign of the SNP, and also the reference and alternative

};

#endif // SNP_H
