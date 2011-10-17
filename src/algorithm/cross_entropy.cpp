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
#include <vector>
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
 * @param[in] screen_output activates output to screen
 * @throws value_error if number of generations is < 1 or elite outside [0,1]
 * 
 * */
cross_entropy::cross_entropy(int gen, double elite, double scale, bool screen_output):base(),
			m_gen(boost::numeric_cast<std::size_t>(gen)),m_elite(elite),
			m_scale(scale),m_screen_output(screen_output) {
	if (gen < 1 || elite < 0 || elite > 1) {
		pagmo_throw(value_error,"number of generation must be > 0 and elite must be in [0,1]");
	}
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
	const problem::base::size_type prob_i_dimension = prob.get_i_dimension(), D = prob.get_dimension(), Dc = D - prob_i_dimension, prob_c_dimension = prob.get_c_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const population::size_type n_elite = boost::numeric_cast<population::size_type>(m_elite * NP);

	//We perform some checks to determine whether the problem/population are suitable for Cross Entropy
	if ( Dc == 0 ) {
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

	if (NP < 2) {
		pagmo_throw(value_error,"for CE at least 2 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}

	if (n_elite < 1) {
		pagmo_throw(value_error,"The population elite contains no individuals ..... maybe increase the elite parameter?");
	}

	using namespace Eigen;
	// We allocate the memory necessary for the multivariate random vector generation
	MatrixXd C(Dc,Dc);					//Covariance Matrix
	MatrixXd U(Dc,Dc);					//Upper Triangular Cholesky
	LLT<MatrixXd> llt(Dc);					//Cholesky Factorization
	VectorXd mu(Dc), tmp(Dc);				//Mean and a temp vector
	std::vector<VectorXd> elite(n_elite,mu);		//Container of the elite chromosomes
	std::vector<population::size_type> elite_idx(n_elite);	//Container of the indexes in pop of the elite
	decision_vector dumb(Dc);

	boost::normal_distribution<double> normal(0,1);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::normal_distribution<double> > normally_distributed_number(m_drng,normal);
	boost::uniform_real<double> uniform(0,1);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > randomly_distributed_number(m_drng,uniform);


	// We start from the champion as the mean
	for (problem::base::size_type i=0;i<Dc;++i){
		mu[i] = pop.champion().x[i];
	}
std::cout << "Starting Mean: " << mu.transpose() << std::endl;
std::cout << "Elite dim: " << n_elite << std::endl<< std::endl;

	// Main loop
	for (std::size_t g = 0; g < m_gen; ++g) {
		// 1 - We extract the elite from this generation (NOTE: we use best_x)
		elite_idx = pop.get_best_idx(n_elite);
		for ( population::size_type i = 0; i<n_elite; i++ ) { 
			for (problem::base::size_type j=0;j<Dc;++j){
				elite[i][j] = pop.get_individual(elite_idx[i]).best_x[j];
			}
		}

		// 2 - We evaluate the Covariance Matrix as least square estimator of the elite (with mean mu)
		tmp = (elite[0] - mu);
		C = tmp*tmp.transpose();
		for ( population::size_type i = 1; i<n_elite; i++ ) { 
			tmp = (elite[i] - mu);
			C = C + tmp*tmp.transpose();
		}
		C = C * (m_scale/n_elite);
std::cout << "C: " << C << std::endl<< std::endl;

		// 3 - We compute the new elite mean
		mu = elite[0];
		for ( population::size_type i = 1; i<n_elite; i++ ) { 
			mu = mu + elite[i];
		}
		mu = mu / n_elite;

std::cout << "Elite Mean: " << mu.transpose() << std::endl;

		// 4 - We sample a new generation
		llt.compute(C);
		U = llt.matrixU();


std::cout << "Upper Triangular: " << U << std::endl  << std::endl ;
		for (population::size_type i = 0; i<NP; i++ ) {
			// 4a - We generate a random vector normally distributed with zero mean and unit variance
			for (problem::base::size_type j=0;j<Dc;++j){
				tmp[j] = normally_distributed_number();
			}
			// 4b - We use Cholesky Triangular form to map the vector into our space

			tmp =  U.transpose()*tmp;
std::cout << "Variation: " << tmp.transpose() << std::endl;
			tmp = tmp+mu;


			// 4c - If generated point is outside the bounds ... fixit!!
			size_t i2 = 0;
			while (i2<Dc) {
				if ((tmp[i2] < lb[i2]) || (tmp[i2] > ub[i2]))
					tmp[i2] = lb[i2] + randomly_distributed_number()*(ub[i2]-lb[i2]);
				++i2;
			}

			for (problem::base::size_type j=0;j<Dc;++j){
				dumb[j] = tmp[j];
			}
			pop.set_x(i,dumb);
		}
	}

}

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
	s << "elite fraction:" << m_elite << ' ';
	s << "covariance scaling:" << m_scale << ' ';
	return s.str();
}



}} //namespaces
