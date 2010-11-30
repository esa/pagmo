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
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>

#include "firefly.h"
#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"




namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] iter number of iterations.
 * @param[in] alpha define the width of the random vector
 * @param[in] beta define the maximum attractiveness
 * @param[in] gamma define the absorption coefficent
 * @throws value_error if number of iterations is negative or alpha, beta and gamma are not in [0,1]
 */
firefly::firefly(int iter, double alpha, double beta, double gamma):base(),m_iter(iter), m_alpha(alpha), m_beta(beta), m_gamma(gamma) {
	if (iter < 0) {
		pagmo_throw(value_error,"number of iterations must be nonnegative");
	}
	if (alpha < 0 || alpha > 1) {
		pagmo_throw(value_error,"alpha should be in [0,1]");
	}
	if (beta < 0 || beta > 1) {
		pagmo_throw(value_error,"beta should be in [0,1]");
	}
	if (gamma < 0 || gamma > 1) {
		pagmo_throw(value_error,"gamma should be in [0,1] intervall");
	}
}

/// Clone method.
base_ptr firefly::clone() const
{
	return base_ptr(new firefly(*this));
}

/// Evolve implementation.
/**
 * Run the Firefly algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void firefly::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type prob_i_dimension = prob.get_i_dimension(), D = prob.get_dimension(), Dc = D - prob_i_dimension, prob_c_dimension = prob.get_c_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = (int) pop.size();

	//We perform some checks to determine wether the problem/population are suitable for Firefly
	if ( Dc == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for Firefly to optimise");
	}

	if ( prob.get_f_dimension() != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and Firefly is not suitable to solve it");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and Firefly is not suitable to solve it");
	}

	if (NP < 2) {
		pagmo_throw(value_error,"for Firefly at least 2 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_iter == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy);	//set of firefly positions
	std::vector<fitness_vector> fit(NP);		//set of firefly positions fitness

	decision_vector::size_type best_firefly = 0;	//best firefly at the begining of an iteration
	fitness_vector best_fitness = fit[0];		//fitness of the best firefly at the begining of an iteration

	double r = 0;	//temp variable to store distances between fireflies
	double b = 0;	//temp variable to store attractiveness of a firefly

	// Copy the fireflies position and their fitness
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;
	}
	
	//Find maximum distance between individuals
	double r_max = 0;
	double r_temp = 0;
	for (population::size_type ii = 0; ii< NP; ++ii) {
		for (population::size_type jj = 0; jj< NP; ++jj) {
			for(problem::base::size_type k=0; k < Dc; ++k) {
				r_temp += (X[ii][k] - X[jj][k])*(X[ii][k]-X[jj][k]);
			}
			r_temp = sqrt(r_temp);
			if (r_temp > r_max) {
				r_max = r_temp;
			}
		}
	}

	double newgamma = m_gamma / r_max;


	// Main Firefly loop
	for (int j = 0; j < m_iter; ++j) {
		
		//Find the best firefly
		best_firefly = 0;
		best_fitness = fit[0];
		for(population::size_type i=1; i < NP; ++i) {
			if(prob.compare_fitness(fit[i], best_fitness)) {
				best_fitness = fit[i];
				best_firefly = i;
			}
		}

		for (population::size_type ii = 0; ii< NP; ++ii) {
			for (population::size_type jj = 0; jj< NP; ++jj) {
				if(prob.compare_fitness(fit[jj], fit[ii])) { //if jj is better than ii
					
					//Calculate distance between X[ii] and X[jj]
					r = 0;
					for(problem::base::size_type k=0; k < Dc; ++k) {
						r += (X[ii][k] - X[jj][k]) * (X[ii][k] - X[jj][k]) ;
					}
					r = sqrt(r);  //distance between X[ii] and X[jj]


					b = m_beta * exp( -1 * newgamma * r*r); //calculate attractiveness

					//Move the firefly ii torwards jj
					for(problem::base::size_type k=0; k < Dc; ++k) {
						X[ii][k] = (1-b) * X[ii][k] + b * X[jj][k] + boost::uniform_real<double>(lb[k]*m_alpha,m_alpha*ub[k])(m_drng);
						
						//check constraints
						if (X[ii][k] < lb[k]) {
							X[ii][k] = lb[k];
						}			
						if (X[ii][k] > ub[k]) {
							X[ii][k] = ub[k];
						}
			
					}
					pop.set_x(ii,X[ii]);
					prob.objfun(fit[ii], X[ii]);
				}
			}

		}
	} // end of main Firefly loop

}

/// Algorithm name
std::string firefly::get_name() const
{
	return "Firefly optimization";
}


/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string firefly::human_readable_extra() const
{
	std::ostringstream s;
	s << "iter:" << m_iter << ' ';
	s << "alpha:" << m_alpha << ' ';
	s << "beta:" << m_beta << ' ';
	s << "gamma:" << m_gamma << ' ';
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::firefly);
