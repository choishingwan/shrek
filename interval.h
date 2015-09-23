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
        Interval(std::string chr, size_t start, size_t end);
        /** Default destructor */
        virtual ~Interval();
        /** Return the chromosome information of this interval */
        inline std::string getChr() const {return m_chr; };
        /** Return the starting loc of this interval */
        inline size_t getStart() const {return m_start; };
        /** Return the last loc of this interval */
        inline size_t getEnd() const { return m_end; };
    protected:
    private:
        std::string m_chr;
        size_t m_start;
        size_t m_end;
};

#endif // INTERVAL_H
