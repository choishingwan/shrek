#ifndef LINKAGE_H
#define LINKAGE_H

class LinkageThread;
#include <deque>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "snp.h"
#include "genotype.h"
#include "linkagethread.h"
#include "processcode.h"


class Linkage
{
	public:
		Linkage(std::vector<Snp*> *snpList, size_t thread);
		virtual ~Linkage();
		ProcessCode Construct(std::deque<Genotype*> &genotype, const size_t &prevResidual, const size_t &blockSize, bool correction, std::vector<size_t> &perfectLd, std::deque<size_t> &snpLoc);
		ProcessCode Initialize(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize);
		ProcessCode Reinitialize(const size_t genotypeSize);
		size_t rows() const;
		size_t cols() const;
		Eigen::MatrixXd block(size_t blockStart, size_t lengthOfBlock);
		Eigen::VectorXd solve(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Eigen::VectorXd *effective);
		double Geteffective() const;
		void Seteffective(double i);
		void Remove(std::vector<size_t> &perfectLd, std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc);
		void print(); //DEBUG
	protected:
	private:
        Eigen::MatrixXd m_linkage;
        std::vector<Snp*> *m_snpList;
        double m_effectiveNumber;
        size_t m_numItem;
        size_t m_thread;
        void triangularThread(const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype, std::vector<size_t> &perfectLd, std::deque<size_t> &snpLoc);
		void rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype, std::vector<size_t> &perfectLd, std::deque<size_t> &snpLoc);
};

#endif // LINKAGE_H
