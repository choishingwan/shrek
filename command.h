#ifndef COMMAND_H
#define COMMAND_H


class Command
{
    public:
        /** Default constructor */
        Command();
        /** Default destructor */
        virtual ~Command();
    protected:
    private:
        size_t m_thread;
        size_t m_blockSize;
        size_t m_sampleSize_;
        size_t m_caseSize_;
        size_t m_controlSize;
        size_t m_cIndex;
        size_t m_tIndex;
        size_t m_bpIndex;
        size_t m_chrIndex;
        size_t m_rsIndex;
        size_t m_sampleSizeIndex;
        double m_prevalence;
        double m_inflation;
        double m_maf;
        bool m_ldCorrection;
        bool m_validate;
        bool m_isPvalue;
        bool m_provideSampleSize;
        bool m_quantitative;
        std::string m_linkageFileName;
        std::string m_outputPrefix;
        std::string m_pValueFileName;
        std::string m_genotypeFilePrefix;
        std::string m_regionList;
};

#endif // COMMAND_H
