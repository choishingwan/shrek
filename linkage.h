#ifndef LINKAGE_H
#define LINKAGE_H

class LinkageThread;
#include <deque>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/IterativeSolvers>
#include <algorithm>
#include <limits>
#include <map>
#include <mutex> //DEBUG
#include <complex>
#include <gsl/gsl_sf_hyperg.h> //debug
#include <gsl/gsl_errno.h> //debug
#include "configure.h"
#include "genotype.h"
#include "linkagethread.h"
#include "processcode.h"
#include "snp.h"

class Linkage
{
	public:
		Linkage();
		//Linkage(size_t thread, std::vector<Snp*> *snpList, std::deque<size_t> *snpLoc);
		virtual ~Linkage();
		ProcessCode Construct(std::deque<Genotype*> &genotype, const size_t &prevResidual, const size_t &blockSize, bool correction);
		ProcessCode Initialize(std::deque<Genotype*> &genotype, const size_t &prevResiduals, const size_t &blockSize);
		ProcessCode Reinitialize(size_t &genotypeSize);
		size_t rows() const;
		size_t cols() const;
		Eigen::MatrixXd block(size_t blockStart, size_t lengthOfBlock);
		Eigen::MatrixXd blockSqrt(size_t blockStart, size_t lengthOfBlock);
		Eigen::MatrixXd varBlock(size_t blockStart, size_t lengthOfBlock);
		Eigen::VectorXd solve(size_t start, size_t length, Eigen::VectorXd *betaEstimate, Eigen::VectorXd *effective);
		Eigen::VectorXd solveChi(size_t start, size_t length, Eigen::VectorXd const* const betaEstimate, Eigen::VectorXd const* const signValue, Eigen::MatrixXd &variance, size_t sampleSize);
		double Geteffective() const;
		static double VarianceR2(double rSq, size_t numSample, size_t predictor);
		static double ExpectedR2(double rSq, size_t numSample, size_t predictor);
        size_t Remove();
        void Update(std::deque<Genotype*> &genotype, std::deque<size_t> &snpLoc);
        void setSnpLoc(std::deque<size_t> *snpLoc);
        void setSnpList(std::vector<Snp* > *snpList);
        void setThread(size_t thread);
		void print();
	protected:
	private:
		static std::mutex mtx;
        Eigen::MatrixXd m_linkage;
        Eigen::MatrixXd m_linkageSqrt;
        Eigen::MatrixXd m_varLinkage;
        size_t m_thread;
        std::vector<size_t> m_perfectLd; //Store the remove index of on matrix level
        std::vector<Snp*> *m_snpList;
        std::deque<size_t> *m_snpLoc;
        void triangularThread(const size_t startBlock, const size_t endBlock, bool correction, std::deque<Genotype*> &genotype);
		void rectangularThread(const size_t start, const size_t width, const size_t height, bool correction, std::deque<Genotype*> &genotype);
};

#endif // LINKAGE_H
