#ifndef LINKAGE_H
#define LINKAGE_H

class LinkageThread;
#include <deque>
#include <Eigen/Dense>
#include <Eigen/SVD>
<<<<<<< HEAD
#include "snp.h"
=======
#include <algorithm>
#include <limits>
#include <map>
#include "configure.h"
>>>>>>> perfectLd
#include "genotype.h"
#include "linkagethread.h"
#include "processcode.h"
#include "snp.h"


class Linkage
{
	public:
<<<<<<< HEAD
		Linkage(std::vector<Snp*> *snpList, size_t thread);
		virtual ~Linkage();
		ProcessCode Construct(std::deque<Genotype*> &genotype, const size_t &prevResidual, const size_t &blockSize, bool correction, std::vector<size_t> &perfectLd, std::deque<size_t> &snpLoc);
		ProcessCode Initialize(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize);
		ProcessCode Reinitialize(const size_t genotypeSize);
=======
		Linkage(size_t thread, std::vector<Snp*> *snpList, std::deque<size_t> *snpLoc);
		virtual ~Linkage();
		ProcessCode Construct(std::deque<Genotype*> &genotype, const size_t &prevResidual, const size_t &blockSize, bool correction);
		ProcessCode Initialize(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize);
		ProcessCode Reinitialize(size_t &genotypeSize);
>>>>>>> perfectLd
		size_t rows() const;
		size_t cols() const;
		Eigen::MatrixXd block(size_t blockStart, size_t lengthOfBlock);
		Eigen::VectorXd solve(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Eigen::VectorXd *effective);
		double Geteffective() const;
<<<<<<< HEAD
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
=======
        size_t Remove();
        void Update(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc);

		void print();
	protected:
	private:
        Eigen::MatrixXd m_linkage;
        size_t m_thread;
        std::vector<size_t> m_perfectLd; //Store the remove index of on matrix level
        std::vector<Snp*> *m_snpList;
        std::deque<size_t> *m_snpLoc;
        void triangularThread(const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype);
		void rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype);
>>>>>>> perfectLd
};

#endif // LINKAGE_H
