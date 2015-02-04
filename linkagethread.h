#ifndef LINKAGETHREAD_H
#define LINKAGETHREAD_H

#include <deque>
#include <mutex>
#include <limits>
#include "linkage.h"
#include "configure.h"
#include "genotype.h"
#include "snp.h"
<<<<<<< HEAD
<<<<<<< HEAD
class LinkageThread
{
	public:
		/** Default constructor */
		LinkageThread(bool correction, const size_t blockEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<Snp*> *snpList, std::vector<size_t> *perfectLd);
		LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<Snp*> *snpList, std::vector<size_t> *perfectLd);
		/** Default destructor */
=======

class LinkageThread
{
	public:
		LinkageThread(bool correction, const size_t blockEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList);
		LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList);
>>>>>>> perfectLd
=======

class LinkageThread
{
	public:
		LinkageThread(bool correction, const size_t blockEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList);
		LinkageThread(bool correction, const size_t snpStart, const size_t snpEnd, const size_t boundStart, const size_t boundEnd, Eigen::MatrixXd *ldMatrix, std::deque<Genotype* > *genotype, std::deque<size_t> *snpLoc, std::vector<size_t> *perfectLd, std::vector<Snp*> *snpList);
>>>>>>> perfectLD
		virtual ~LinkageThread();

		void Addstart(size_t i);
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
        std::deque<size_t> *m_snpLoc;
<<<<<<< HEAD
<<<<<<< HEAD
        std::vector<Snp*> *m_snpList;
        std::vector<size_t> *m_perfectLd;
        std::vector<size_t> m_startLoc;
        void triangularProcess();
        void rectangularProcess();
=======
=======
>>>>>>> perfectLD
        std::vector<size_t> m_startLoc;
        std::vector<size_t> *m_perfectLd;
        std::vector<Snp*> *m_snpList;
		void triangularProcess();
        void rectangularProcess();
        static std::mutex mtx;
<<<<<<< HEAD
>>>>>>> perfectLd
=======
>>>>>>> perfectLD
};

#endif // LINKAGETHREAD_H
