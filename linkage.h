#ifndef LINKAGE_H
#define LINKAGE_H

#include <deque>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "genotype.h"

class Linkage
{
	public:
		Linkage();
		virtual ~Linkage();
		void Construct(std::deque<Genotype*> genotype);
	protected:
	private:
        Eigen::MatrixXd m_linkage;
};

#endif // LINKAGE_H
