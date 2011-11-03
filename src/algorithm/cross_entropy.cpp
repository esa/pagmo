/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include <string>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <math.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/numeric/conversion/cast.hpp>


#include "cross_entropy.h"
#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "../Eigen/Dense"




namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations
 * @param[in] elite the fraction of samples to be considered elite
 * @param[in] scale multiplication coefficient for the generated points
 * @param[in] variant algoritmic variant to use (one of [1,2])
		       1. 'Canonical' - Covariance Matrix is evaluated as sum (x_(i+1)-mu_i)^T (x_(i+1)-mu_i)
		       2. 'Dario's' - Covariance Matrix is evaluated as   sum (x_(i+1)-mu_i^T)^T (x_(i+1)-mu_i^T)
 * @param[in] screen_output activates output to screen
 * @throws value_error if number of generations is < 1 or elite outside [0,1]
 * 
 * */
cross_entropy::cross_entropy(int gen, double cc, double cs, double c1, double cmu, double sigma0, double ftol, double xtol):base(),
			m_gen(boost::numeric_cast<std::size_t>(gen)), m_cc(cc), m_cs(cs), m_c1(c1), m_cmu(cmu), m_sigma(sigma0),
			m_ftol(ftol), m_xtol(xtol), m_screen_output(false) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if ( ((cc < 0) || (cc > 0)) && !(cc==-1) ){
		pagmo_throw(value_error,"cc needs to be in [0,1] or -1 for auto value");
	}
	if ( ((cs < 0) || (cs > 0)) && !(cs==-1) ){
		pagmo_throw(value_error,"cs needs to be in [0,1] or -1 for auto value");
	}
	if ( ((c1 < 0) || (c1 > 0)) && !(c1==-1) ){
		pagmo_throw(value_error,"c1 needs to be in [0,1] or -1 for auto value");
	}
	if ( ((cmu < 0) || (cmu > 0)) && !(cmu==-1) ){
		pagmo_throw(value_error,"cmu needs to be in [0,1] or -1 for auto value");
	}

	//Initialize the algorithm memory
	m_mean = Eigen::VectorXd(1);
	m_B = Eigen::MatrixXd(1,1);
	m_D = Eigen::VectorXd(1);
	m_C = Eigen::MatrixXd(1,1);
	m_invsqrtC = Eigen::MatrixXd(1,1);
	m_pc = Eigen::VectorXd(1);
	m_ps = Eigen::VectorXd(1);
	m_counteval = 0;
	m_eigeneval = 0;

}
/// Clone method.
base_ptr cross_entropy::clone() const
{
	return base_ptr(new cross_entropy(*this));
}

/// Evolve implementation.
/**
 * Run Cross Entropy
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void cross_entropy::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type prob_i_dimension = prob.get_i_dimension(), dim = prob.get_dimension(), N = dim - prob_i_dimension, prob_c_dimension = prob.get_c_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type lam = pop.size();
	const population::size_type mu = boost::numeric_cast<population::size_type>(lam/2);

	//We perform some checks to determine whether the problem/population are suitable for Cross Entropy
	if ( N == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for CE to optimise");
	}

	if ( prob.get_f_dimension() != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and CE is not suitable to solve it");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and CE is not suitable to solve it");
	}

	if ( prob_i_dimension != 0 ) {
		pagmo_throw(value_error,"The problem has an integer part and CE is not suitable to solve it");
	}

	if (lam < 5) {
		pagmo_throw(value_error,"for CE at least 5 individuals in the population are required");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}

	using namespace Eigen;

	// Initializing the random number generators
	boost::normal_distribution<double> normal(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::normal_distribution<double> > normally_distributed_number(m_drng,normal);
	boost::uniform_real<double> uniform(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > randomly_distributed_number(m_drng,uniform);

	// Setting coefficients for Selection
	VectorXd weights(mu);
	for (int i = 0; i < weights.rows(); ++i){
		weights[i] = log(mu+0.5) - log(i+1);
	}
	weights /= weights.sum();						// weights for weighted recombination
	double mueff = 1.0 / (weights.transpose()*weights);			// variance-effectiveness of sum w_i x_i

	// Setting coefficients for Adaptation automatically or to user defined data
	double cc(m_cc), cs(m_cs), c1(m_c1), cmu(m_cmu);
	if (m_cc == -1) {
		cc = (4 + mueff/N) / (N+4 + 2*mueff/N);			// t-const for cumulation for C
	}
	if (m_cs == -1) {
		cs = (mueff+2) / (N+mueff+5);				// t-const for cumulation for sigma control
	}
	if (m_c1 == -1) {
		c1 = 2 / ((N+1.3)*(N+1.3)+mueff);			// learning rate for rank-one update of C
	}
	if (m_cmu == -1) {
		cmu = 2 * (mueff-2+1/mueff) / ((N+2)*(N+2)+mueff);	// and for rank-mu update
	}
	
	double damps = 1 + 2*std::max(0.0, sqrt((mueff-1)/(N+1))-1) + cs;	// damping for sigma
	double chiN = pow(N,0.5)*(1-1/(4*N)+1/(21*N*N));			// expectation of ||N(0,I)|| == norm(randn(N,1))

	// Initializing and allocating (here one could use mutable data member to avoid redefinition of non const data)
	VectorXd mean(m_mean);
	VectorXd meanold(N);
	VectorXd variation(m_variation);
	std::vector<VectorXd> newpop(m_newpop);
	MatrixXd B(m_B);
	MatrixXd D(m_D);
	MatrixXd Dinv(N,N);
	MatrixXd C(m_C);
	MatrixXd Cold(N,N);
	MatrixXd invsqrtC(m_invsqrtC);
	VectorXd pc(m_pc);
	VectorXd ps(m_ps);
	int counteval(m_counteval);
	int eigeneval(m_eigeneval);
	double sigma(m_sigma);

	VectorXd tmp = VectorXd::Zero(N);
	std::vector<VectorXd> elite(mu,tmp);

	decision_vector dumb(N,0);

	// If the algorithm is called for the first time on this problem dimension / pop size we erease the memory
	if ( !((m_newpop.size() == lam) && (m_newpop[0].rows() == N)) ) {
		for (problem::base::size_type i=0;i<N;++i){
			mean(i) = pop.champion().x[i];
		}

		newpop = std::vector<VectorXd>(lam,tmp);
		variation = VectorXd(N);
		B = MatrixXd::Identity(N,N);					//B defines the coordinate system
		D = MatrixXd::Identity(N,N);					//diagonal D defines the scaling
		C = MatrixXd::Identity(N,N);					//covariance matrix C
		invsqrtC = MatrixXd::Identity(N,N);					//inverse of sqrt(C)
		pc = VectorXd::Zero(N);
		ps = VectorXd::Zero(N);
		counteval = 0;
		eigeneval = 0;
	}
	
	// ----------------------------------------------//
	// HERE WE START THE REAL ALGORITHM              //
	// ----------------------------------------------//
	
	SelfAdjointEigenSolver<MatrixXd> es(N);
	for (std::size_t g = 0; g < m_gen; ++g) {
		// 1 - We generate and evaluate lam new individuals
		for (population::size_type i = 0; i<lam; ++i ) {
			// 1a - we create a randomly normal distributed vector
			for (problem::base::size_type j=0;j<N;++j){
				tmp(j) = normally_distributed_number();
			}
			// 1b - and store its transformed value in the newpop
			newpop[i] = mean + sigma * B * D * tmp;
		}
		// 1c - we fix the bounds and reinsert
		for (population::size_type i = 0; i<lam; ++i ) {
			for (decision_vector::size_type j = 0; j<N; ++i ) {
				if ((newpop[i](j) < lb[j]) || (newpop[i](j) > ub[j]) ) {
					newpop[i](j) = lb[j] + randomly_distributed_number() * (ub[j] - lb[j]);
				}
			}
		}
		for (population::size_type i = 0; i<lam; ++i ) {
			for (decision_vector::size_type j = 0; j<N; ++i ) {
				dumb[j] = newpop[i](j);
			}
			int idx = pop.get_worst_idx();
			pop.set_x(idx,dumb);
		}
		counteval += lam;

		// 2 - We extract the elite from this generation
		std::vector<population::size_type> best_idx = pop.get_best_idx(mu);
		for (population::size_type i = 0; i<mu; ++i ) {
			for (decision_vector::size_type j = 0; j<N; ++i ) {
				elite[i](j) = pop.get_individual(best_idx[i]).best_x[j];
			}
		}

		// 3 - Compute the new elite mean storing the old one
		meanold=mean;
		mean = elite[0]*weights[0];
		for (population::size_type i = 0; i<mu; ++i ) {
			mean += elite[i]*weights[i];
		}

		// 4 - Update evolution paths
		ps = (1 - cs) * ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (mean-meanold) / sigma;
		int hsig = (ps.squaredNorm() / N / (1-pow((1-cs),(2*counteval/lam))) ) < (2 + 4/(N+1));
		pc = (1-cc) * pc + hsig * sqrt(cc*(2-cc)*mueff) * (mean-meanold) / sigma;

		// 5 - Adapt Covariance Matrix
		Cold = C;
		C = (elite[0]-meanold)*(elite[0]-meanold).transpose()*weights[0];
		for (population::size_type i = 0; i<mu; ++i ) {
			C += (elite[i]-meanold)*(elite[i]-meanold).transpose()*weights[i];
		}
		C /= sigma*sigma;
		C = (1-c1-cmu)*Cold + cmu*C + c1 * ((pc * pc.transpose()) + (1-hsig) * cc*(2-cc) * Cold);

		//6 - Adapt sigma
		sigma *= exp( (cs/damps)*(ps.norm()/chiN - 1));

		//7 - Perform eigen-decomposition of C
		if ( (counteval - eigeneval) > (lam/(c1+cmu)/N/10) ) {		//achieve O(N^2)
			eigeneval = counteval;
			//C = (C+C.T)/2						//enforce symmetry
			es.solve(C);						//eigen decomposition
			B = es.eigenvectors();
			D = es.eigenvalues().asDiagonal();
			for (decision_vector::size_type j = 0; j<N; ++i ) {
				D(j,j) = sqrt(D(j,j));				//D contains standard deviations now
			}
			for (decision_vector::size_type j = 0; j<N; ++i ) {
				Dinv(j,j) = 1 / D(j,j);
			}
			invsqrtC = B*Dinv*B.transpose();
		}
	}
	


	
		
	// 6 - We print on screen if required
	if (m_screen_output) {
		if (1)
			std::cout << std::endl << std::left << std::setw(20) <<"Gen." << std::setw(20) << "Champion " << 
				std::setw(20) << "Best " << std::setw(20) << "Worst" << std::setw(20) << "Variation" << std::endl; 
				
		std::cout << std::left << std::setprecision(14) << std::setw(20) << 
			1 << std::setw(20)<< pop.champion().f[0] << std::setw(20) << 
			pop.get_individual(pop.get_best_idx()).best_f[0] << std::setw(20) << 
			pop.get_individual(pop.get_worst_idx()).best_f[0] << std::setw(20) << 32.22 << std::endl;
		
	}

}

/// Sets screen output
/**
 * Sets CES screen output at a default level
 *
 * @param[in] p true or false
 */
void cross_entropy::set_screen_output(const bool p) {m_screen_output = p;}

/// Gets screen output
/**
 * Gets CES screen output level
 *
 * @param[out] boolean 
 */
bool cross_entropy::get_screen_output() const {return m_screen_output;}

/// Algorithm name
std::string cross_entropy::get_name() const
{
	return "Cross Entropy Study";
}


/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string cross_entropy::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	return s.str();
}



}} //namespaces
