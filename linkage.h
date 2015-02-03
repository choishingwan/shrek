#ifndef LINKAGE_H
#define LINKAGE_H

class LinkageThread;
#include <deque>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <algorithm>
#include <limits>
#include <map>
#include "configure.h"
#include "genotype.h"
#include "linkagethread.h"
#include "processcode.h"
#include "snp.h"


class Linkage
{
	public:
		Linkage(size_t thread);
		virtual ~Linkage();
		ProcessCode Construct(std::deque<Genotype*> &genotype, const size_t &prevResidual, const size_t &blockSize, bool correction);
		size_t rows() const;
		size_t cols() const;
		Eigen::MatrixXd block(size_t blockStart, size_t lengthOfBlock);
		Eigen::VectorXd solve(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Eigen::VectorXd *effective);
		double Geteffective() const;
        size_t Remove();
        void Update(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc, std::vector<Snp*> *snpList);
		void Seteffective(double i);
		void print();
	protected:
	private:
        Eigen::MatrixXd m_linkage;
        double m_effectiveNumber;
        size_t m_numItem;
        size_t m_thread;
        std::vector<size_t> m_perfectLd;
        void triangularThread(const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype);
		void rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype);
};

#endif // LINKAGE_H
