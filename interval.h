#ifndef INTERVAL_H
#define INTERVAL_H

#include <string>

/** \class Interval
 *  \brief Storage class for coordinates
 */
class Interval
{
    public:
        /** Default constructor */
        Interval(std::string chr, size_t start, size_t end):m_chr(chr), m_start(start), m_end(end){};
        /** Default destructor */
        virtual ~Interval(){ };
        /** Return the chromosome information of this interval */
        std::string getChr() const {return m_chr; };
        /** Return the starting loc of this interval */
        size_t getStart() const {return m_start; };
        /** Return the last loc of this interval */
        size_t getEnd() const { return m_end; };
        void setEnd(size_t i ){ m_end = i;};
    protected:
    private:
        std::string m_chr="";
        size_t m_start=0;
        size_t m_end=0;
};

#endif // INTERVAL_H
