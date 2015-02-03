#ifndef LINKAGETHREAD_H
#define LINKAGETHREAD_H

#include <deque>
#include <mutex>
#include <limits>
#include "linkage.h"
#include "configure.h"
#include "genotype.h"

class LinkageThread
{
	public:
		/** Default constructor */
		LinkageThread(bool correction, const size_t blockEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::vector<size_t> *perfectLd);
		LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::vector<size_t> *perfectLd);
		/** Default destructor */
		virtual ~LinkageThread();

		void Addstart(size_t i);
		bool Getcorrection() const;
		size_t GetsnpStart() const;
		size_t GetsnpEnd() const;
		size_t GetboundStart() const;
		size_t GetboundEnd() const;
		size_t GetstartLoc(size_t i) const;
		size_t GetsizeOfStart() const;
		Eigen::MatrixXd *Getld();
		std::deque<Genotype* > *Getgenotype();
        static void *triangularProcess(void *in);
        static void *rectangularProcess(void *in);
	protected:
	private:
        bool m_correction;
        size_t m_snpStart;
        size_t m_snpEnd;
        size_t m_boundStart;
        size_t m_boundEnd;
        Eigen::MatrixXd *m_ldMatrix;
        std::deque<Genotype*> *m_genotype;
        std::vector<size_t> m_startLoc;
        std::vector<size_t> *m_perfectLd;
		void triangularProcess();
        void rectangularProcess();
        static std::mutex mtx;
};

#endif // LINKAGETHREAD_H
