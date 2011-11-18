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
#include <cmath>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/numeric/conversion/cast.hpp>


#include "cmaes.h"
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
 * @param[in] cc time constant for C cumulation (in [0,1]) if -1 automatic values are set
 * @param[in] cs time constant for sigma cumulation (in [0,1]) if -1 automatic values are set
 * @param[in] c1 learning rate for rank-1 update (in [0,1]) if -1 automatic values are set
 * @param[in] cmu learning rate for rank-mu update (in [0,1]) if -1 automatic values are set
 * @param[in] sigma0 starting step (std)
 * @param[in] ftol stopping criteria on the x tolerance
 * @param[in] xtol stopping criteria on the f tolerance
 * @param[in] restart when the algorithm looses its memory of the parameter adaptation (C, p etc ....) at each call
 * @throws value_error if cc,cs,c1,cmu are not in [0,1] or not -1
 * 
 * */
cmaes::cmaes(int gen, double cc, double cs, double c1, double cmu, double sigma0, double ftol, double xtol, bool restart):
		base(), m_gen(boost::numeric_cast<std::size_t>(gen)), m_cc(cc), m_cs(cs), m_c1(c1), 
		m_cmu(cmu), m_sigma(sigma0), m_ftol(ftol), m_xtol(xtol), m_restart(restart) {
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
	m_mean = Eigen::VectorXd::Zero(1);
	m_variation = Eigen::VectorXd::Zero(1);
	m_newpop = std::vector<Eigen::VectorXd>();
	m_B = Eigen::MatrixXd::Identity(1,1);
	m_D = Eigen::MatrixXd::Identity(1,1);
	m_C = Eigen::MatrixXd::Identity(1,1);
	m_invsqrtC = Eigen::MatrixXd::Identity(1,1);
	m_pc = Eigen::VectorXd::Zero(1);
	m_ps = Eigen::VectorXd::Zero(1);
	m_counteval = 0;
	m_eigeneval = 0;

}
/// Clone method.
base_ptr cmaes::clone() const
{
	return base_ptr(new cmaes(*this));
}

/// Evolve implementation.
/**
 * Run Cross Entropy
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void cmaes::evolve(population &pop) const
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
		weights(i) = std::log(mu+0.5) - std::log(i+1);
	}
	weights /= weights.sum();					// weights for weighted recombination
	double mueff = 1.0 / (weights.transpose()*weights);		// variance-effectiveness of sum w_i x_i

	// Setting coefficients for Adaptation automatically or to user defined data
	double cc(m_cc), cs(m_cs), c1(m_c1), cmu(m_cmu);
	if (cc == -1) {
		cc = (4 + mueff/N) / (N+4 + 2*mueff/N);			// t-const for cumulation for C
	}
	if (cs == -1) {
		cs = (mueff+2) / (N+mueff+5);				// t-const for cumulation for sigma control
	}
	if (c1 == -1) {
		c1 = 2.0 / ((N+1.3)*(N+1.3)+mueff);			// learning rate for rank-one update of C
	}
	if (cmu == -1) {
		cmu = 2.0 * (mueff-2+1/mueff) / ((N+2)*(N+2)+mueff);	// and for rank-mu update
	}
	
	double damps = 1 + 2*std::max(0.0, std::sqrt((mueff-1)/(N+1))-1) + cs;	// damping for sigma
	double chiN = std::sqrt(N) * (1-1.0/(4*N)+1.0/(21*N*N));		// expectation of ||N(0,I)|| == norm(randn(N,1))

	// Initializing and allocating (here one could use mutable data member to avoid redefinition of non const data)

	// Algorithm's Memory. This allows the algorithm to start from its last "state"
	VectorXd mean(m_mean);
	VectorXd variation(m_variation);
	std::vector<VectorXd> newpop(m_newpop);
	MatrixXd B(m_B);
	MatrixXd D(m_D);
	MatrixXd C(m_C);
	MatrixXd invsqrtC(m_invsqrtC);
	VectorXd pc(m_pc);
	VectorXd ps(m_ps);
	int counteval(m_counteval);
	int eigeneval(m_eigeneval);
	double sigma(m_sigma);
	double var_norm = 0;

	// Some buffers
	VectorXd meanold = VectorXd::Zero(N);
	MatrixXd Dinv = MatrixXd::Identity(N,N);
	MatrixXd Cold = MatrixXd::Identity(N,N);
	VectorXd tmp = VectorXd::Zero(N);
	std::vector<VectorXd> elite(mu,tmp);
	decision_vector dumb(N,0);

	// If the algorithm is called for the first time on this problem dimension / pop size or if m_fresh_start is true we erease the memory of past calls
	if ( (m_newpop.size() != lam) || ((unsigned int)(m_newpop[0].rows() ) != N) || (m_restart==true) ) {
		mean.resize(N);
		for (problem::base::size_type i=0;i<N;++i){
			mean(i) = pop.champion().x[i];
		}
		newpop = std::vector<VectorXd>(lam,tmp);
		variation.resize(N);
		B.resize(N,N); B = MatrixXd::Identity(N,N);			//B defines the coordinate system
		D.resize(N,N); D = MatrixXd::Identity(N,N);			//diagonal D defines the scaling
		C.resize(N,N); C = MatrixXd::Identity(N,N);			//covariance matrix C
		invsqrtC.resize(N,N); invsqrtC = MatrixXd::Identity(N,N);	//inverse of sqrt(C)
		pc.resize(N); pc = VectorXd::Zero(N);
		ps.resize(N); ps = VectorXd::Zero(N);
		counteval = 0;
		eigeneval = 0;
	}
	
	// ----------------------------------------------//
	// HERE WE START THE REAL ALGORITHM              //
	// ----------------------------------------------//

	if (m_screen_output) {
		std::cout << "CMAES 4 PaGMO: " << std::endl;
		std::cout << "mu: " << mu
			<< " - lambda: " << lam
			<< " - mueff: " << mueff
			<< " - N: " << N << std::endl;

		std::cout << "cc: " << cc
			<< " - cs: " << cs
			<< " - c1: " << c1
			<< " - cmu: " << cmu
			<< " - sigma: " << sigma
			<< " - damps: " << damps
			<< " - chiN: " << chiN << std::endl;
	}
	
	SelfAdjointEigenSolver<MatrixXd> es(N);
	for (std::size_t g = 0; g < m_gen; ++g) {
		// 1 - We generate and evaluate lam new individuals

		for (population::size_type i = 0; i<lam; ++i ) {
			// 1a - we create a randomly normal distributed vector
			for (problem::base::size_type j=0; j<N; ++j){
				tmp(j) = normally_distributed_number();
			}
			// 1b - and store its transformed value in the newpop
			newpop[i] = mean + sigma * B * D * tmp;
		}
		var_norm = (sigma * B * D * tmp).norm();

		// 1c - we fix the bounds and reinsert
		for (population::size_type i = 0; i<lam; ++i ) {
			for (decision_vector::size_type j = 0; j<N; ++j ) {
				if ( (newpop[i](j) < lb[j]) || (newpop[i](j) > ub[j]) ) {
					newpop[i](j) = lb[j] + randomly_distributed_number() * (ub[j] - lb[j]);
				}
			}
		}

		for (population::size_type i = 0; i<lam; ++i ) {
			for (decision_vector::size_type j = 0; j<N; ++j ) {
				dumb[j] = newpop[i](j);
			}
			int idx = pop.get_worst_idx();
			pop.set_x(idx,dumb);
		}
		counteval += lam;

		// 2 - We extract the elite from this generation
		std::vector<population::size_type> best_idx = pop.get_best_idx(mu);
		for (population::size_type i = 0; i<mu; ++i ) {
			for (decision_vector::size_type j = 0; j<N; ++j ) {
				elite[i](j) = pop.get_individual(best_idx[i]).best_x[j];
			}
		}

		// 3 - Compute the new elite mean storing the old one
		meanold=mean;
		mean = elite[0]*weights(0);
		for (population::size_type i = 1; i<mu; ++i ) {
			mean += elite[i]*weights(i);
		}

		// 4 - Update evolution paths
		ps = (1 - cs) * ps + std::sqrt(cs*(2-cs)*mueff) * invsqrtC * (mean-meanold) / sigma;
		double hsig = 0;
		hsig = (ps.squaredNorm() / N / (1-std::pow((1-cs),(2.0*counteval/lam))) ) < (2.0 + 4/(N+1));
		pc = (1-cc) * pc + hsig * std::sqrt(cc*(2-cc)*mueff) * (mean-meanold) / sigma;

		// 5 - Adapt Covariance Matrix
		Cold = C;
		C = (elite[0]-meanold)*(elite[0]-meanold).transpose()*weights(0);
		for (population::size_type i = 1; i<mu; ++i ) {
			C += (elite[i]-meanold)*(elite[i]-meanold).transpose()*weights(i);
		}
		C /= sigma*sigma;
		C = (1-c1-cmu) * Cold +
			cmu * C +
			c1 * ((pc * pc.transpose()) + (1-hsig) * cc * (2-cc) * Cold);

		//6 - Adapt sigma
		sigma *= std::exp( (cs/damps)*(ps.norm()/chiN - 1) );
		if ( std::isnan(sigma) || std::isinf(sigma) || std::isinf(var_norm) || std::isnan(var_norm) ) {
			std::cout << "eigen: " << es.info() << std::endl;
			std::cout << "B: " << B << std::endl;
			std::cout << "D: " << D << std::endl;
			std::cout << "Dinv: " << D << std::endl;
			std::cout << "invsqrtC: " << invsqrtC << std::endl;
			pagmo_throw(value_error,"NaN!!!!!");
		}

		//7 - Perform eigen-decomposition of C
		if ( (counteval - eigeneval) > (lam/(c1+cmu)/N/10) ) {		//achieve O(N^2)
			eigeneval = counteval;
			C = (C+C.transpose())/2;				//enforce symmetry
			es.compute(C);						//eigen decomposition
			if (es.info()==Success) {
				B = es.eigenvectors();
				D = es.eigenvalues().asDiagonal();
				for (decision_vector::size_type j = 0; j<N; ++j ) {
					D(j,j) = std::sqrt( std::max(1e-20,D(j,j)) );				//D contains standard deviations now
				}
				for (decision_vector::size_type j = 0; j<N; ++j ) {
					Dinv(j,j) = 1.0 / D(j,j);
				}
				invsqrtC = B*Dinv*B.transpose();
			} //if eigendecomposition fails just skip it and keep pevious succesful one.
		}

		//8 - We print on screen if required
		if (m_screen_output) {
			if (!(g%20)) {
				std::cout << std::endl << std::left << std::setw(20) << 
				"Gen." << std::setw(20) << 
				"Champion " << std::setw(20) << 
				"Highest " << std::setw(20) << 
				"Lowest" << std::setw(20) << 
				"Variation" << std::setw(20) << 
				"Step" << std::endl; 
			}
				
			std::cout << std::left << std::setprecision(14) << std::setw(20) << 
				g << std::setw(20) << 
				pop.champion().f[0] << std::setw(20) << 
				pop.get_individual(pop.get_best_idx()).best_f[0] << std::setw(20) << 
				pop.get_individual(pop.get_worst_idx()).best_f[0] << std::setw(20) << 
				var_norm << std::setw(20) <<
				sigma << std::endl;
		}

		//9 - Check the exit conditions (every 40 generations)
		if (g%40) {
			if  ( (sigma*B*D*tmp).norm() < m_xtol ) {
				if (m_screen_output) { 
					std::cout << "Exit condition -- xtol < " <<  m_xtol << std::endl;
				}
				return;
			}

			double mah = std::fabs(pop.get_individual(pop.get_worst_idx()).best_f[0] - pop.get_individual(pop.get_best_idx()).best_f[0]);

			if (mah < m_ftol) {
				if (m_screen_output) {
					std::cout << "Exit condition -- ftol < " <<  m_ftol << std::endl;
				}
				return;
			}
		}

	// Update algorithm memory
	if (!m_restart) {
		m_mean = mean;
		m_variation = variation;
		m_newpop = newpop;
		m_B = B;
		m_D = D;
		m_C = C;
		m_invsqrtC = invsqrtC;
		m_pc = pc;
		m_ps = ps;
		m_counteval = counteval;
		m_eigeneval = eigeneval;
		m_sigma = sigma;
	}
		
	} // end loop on g
}


void cmaes::set_gen(const int gen) {m_gen = gen;}
int cmaes::get_gen() const {return m_gen;}

void cmaes::set_cc(const double cc) {m_cc = cc;}
double cmaes::get_cc() const {return m_cc;}

void cmaes::set_cs(const double cs) {m_cs = cs;}
double cmaes::get_cs() const {return m_cs;}

void cmaes::set_c1(const double c1) {m_c1 = c1;}
double cmaes::get_c1() const {return m_c1;}

void cmaes::set_cmu(const double cmu) {m_cmu = cmu;}
double cmaes::get_cmu() const {return m_cmu;}

void cmaes::set_sigma(const double sigma) {m_sigma = sigma;}
double cmaes::get_sigma() const {return m_sigma;}

void cmaes::set_ftol(const double ftol) {m_ftol = ftol;}
double cmaes::get_ftol() const {return m_ftol;}

void cmaes::set_xtol(const double xtol) {m_xtol = xtol;}
double cmaes::get_xtol() const {return m_xtol;}

/// Algorithm name
std::string cmaes::get_name() const
{
	return "CMAES";
}


/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string cmaes::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' '
	  << "cc:" << m_cc << ' '
	  << "cs:" << m_cs << ' '
	  << "c1:" << m_c1 << ' '
	  << "cmu:" << m_cmu << ' '
	  << "sigma0:" << m_sigma << ' '
	  << "ftol:" << m_ftol << ' '
	  << "xtol:" << m_xtol << ' ' 
	  << "restart:" << m_restart;
	return s.str();
}



}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cmaes);
