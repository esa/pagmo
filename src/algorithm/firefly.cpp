/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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
 * @param[in] gen number of iterations.
 * @param[in] alpha define the width of the random vector
 * @param[in] beta define the maximum attractiveness
 * @param[in] gamma define the absorption coefficent
 * @throws value_error if number of iterations is negative or alpha, beta and gamma are not in [0,1]
 */
firefly::firefly(int gen, double alpha, double beta, double gamma):base(),m_iter(gen), m_alpha(alpha), m_beta(beta), m_gamma(gamma) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of iterations must be nonnegative");
	}
	if (alpha < 0 || alpha > 1) {
		pagmo_throw(value_error,"alpha should be in [0,1]");
	}
	if (beta < 0 || beta > 1) {
		pagmo_throw(value_error,"beta should be in [0,1]");
	}
	if (gamma < 0 || gamma > 1) {
		pagmo_throw(value_error,"gamma should be in [0,1] interval");
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

	//We perform some checks to determine whether the problem/population are suitable for Firefly
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
	std::vector<decision_vector> X(NP,dummy);	//set of firefly positions
	std::vector<decision_vector> X0(NP,dummy);	//set of firefly positions kept to calculate velocity
	std::vector<fitness_vector> fit(NP);		//set of firefly positions fitness

	// Copy the fireflies position and their fitness
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;
		X0[i]	=	pop.get_individual(i).cur_x;
	}
	
	double gamma_nominal_distance = 16.0; // factor comes from scaling r_sqrd by r_max_sqrd in attractiveness calculation (applies a nominal distance)
        double newgamma = gamma_nominal_distance * m_gamma;

	// Main Firefly loop
	for (int j = 0; j < m_iter; ++j) {

		//Find maximum distance between individuals
		double r_max_sqrd = 0;
		for (population::size_type ii = 0; ii< NP; ++ii) {
			for (population::size_type jj = ii+1; jj< NP; ++jj) {
				double r_temp_sqrd = 0;
				for(problem::base::size_type k=0; k < Dc; ++k) {
					r_temp_sqrd += (X[ii][k] - X[jj][k])*(X[ii][k]-X[jj][k]);
				}
				if (r_temp_sqrd > r_max_sqrd) {
					r_max_sqrd = r_temp_sqrd;
				}
			}
		}

		bool moveIItoJJ;
		double r_sqrd;    //temp variable to store distances squared between fireflies
		double b;    //temp variable to store attractiveness of a firefly
		decision_vector X_start;
		fitness_vector test_fit = fit[0];
		for (population::size_type ii = 0; ii< NP; ++ii) {
			for (population::size_type jj = 0; jj< NP; ++jj) {
				moveIItoJJ = prob.compare_fitness(fit[jj], fit[ii]);    //if jj is better than ii
				if(moveIItoJJ) { 

					//Calculate distance between X[ii] and X[jj]
					r_sqrd = 0;
					for(problem::base::size_type k=0; k < Dc; ++k) {
						r_sqrd += (X[ii][k] - X[jj][k]) * (X[ii][k] - X[jj][k]) ;
					}

					b = m_beta * exp( -1 * newgamma * sqrt(r_sqrd/r_max_sqrd)); //calculate attractiveness

                                        //Move the firefly ii torwards jj
					for(problem::base::size_type k=0; k < Dc; ++k) {
						X[ii][k] = (1-b) * X[ii][k] + b * X[jj][k];
					}
				}
				else {
					X_start = X[ii];
				}

				// always apply random walk and check bounds
				for(problem::base::size_type k = 0; k< Dc; ++k) {
					X[ii][k] += boost::uniform_real<double>(-m_alpha, m_alpha)(m_drng) * (ub[k] - lb[k]);
					
					//check constraints
					if (X[ii][k] < lb[k]) {
						X[ii][k] = lb[k];
					}			
					else if (X[ii][k] > ub[k]) {
						X[ii][k] = ub[k];
					}
			
				}

				prob.objfun(test_fit, X[ii]);
				if(moveIItoJJ || prob.compare_fitness(test_fit, fit[ii])) { // only if moving ii towards jj or if new location has better fitness, update population and fitness
					pop.set_x(ii, X[ii]);
                                        fit[ii] = test_fit;
				}
				else {
					X[ii] = X_start;
				}
			}

		}
	} // end of main Firefly loop
	for (population::size_type i = 0; i< NP; ++i) {
		std::transform(X[i].begin(), X[i].end(), X0[i].begin(), dummy.begin(), std::minus<double>()); // dummy is now velocity for i-th individual
		pop.set_v(i, dummy);
	}

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

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::firefly)
