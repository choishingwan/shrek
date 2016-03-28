#ifndef INTERVAL_H
#define INTERVAL_H

#include <string>

/** \class Interval
 *  \brief Storage class for coordinates
 */
class Interval
{
    public:
        Interval(std::string chr, size_t start, size_t end):m_chr(chr),m_start(start),m_end(end){};
        virtual ~Interval(){};
        std::string getChr() const {return m_chr; };
        size_t getStart() const {return m_start; };
        size_t getEnd() const { return m_end; };
        void setEnd(size_t i ){ m_end = i;};
        bool operator < (const Interval &a) const{
            if(m_start==a.m_start) return m_end < a.m_end;
            return m_start <a.m_start;
        }
    protected:
    private:
        std::string m_chr;
        size_t m_start;
        size_t m_end;
};

#endif // INTERVAL_H
