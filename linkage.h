#ifndef LINKAGE_H
#define LINKAGE_H

class LinkageThread;
#include <deque>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "genotype.h"
#include "linkagethread.h"
#include "processcode.h"


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
		double Getbug() const; //DEBUG
		void Setbug(double i) { m_bug+= i; } //DEBUG
	protected:
	private:
        Eigen::MatrixXd m_linkage;
        double m_effectiveNumber;
        size_t m_numItem;
        size_t m_thread;
        double m_bug; // DEBUG
        void triangularThread(const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype);
		void rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype);
};

#endif // LINKAGE_H
