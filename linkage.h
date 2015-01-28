#ifndef LINKAGE_H
#define LINKAGE_H

class LinkageThread;
#include <deque>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "genotype.h"
#include "linkagethread.h"

class Linkage
{
	public:
		Linkage(size_t thread);
		virtual ~Linkage();
		void Construct(std::deque<Genotype*> &genotype, const size_t &prevResidual, const size_t &blockSize, bool correction);
	protected:
	private:
        Eigen::MatrixXd m_linkage;
        double m_effectiveNumber;
        size_t m_numItem;
        size_t m_thread;
        void triangularThread(const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype);
		void rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype);
};

#endif // LINKAGE_H
