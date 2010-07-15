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
#include <boost/numeric/conversion/cast.hpp>

#include "aco.h"
#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../problem/base_aco.h"
#include "../types.h"




namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] iter number of iterations.
 * @param[in] rho evaporation rate
 * @throws value_error if number of iterations is negativ
 */
aco::aco(int iter, double rho):base(),m_iter(iter),m_rho(rho) {
	if (iter < 0) {
		pagmo_throw(value_error,"number of iterations must be nonnegative");
	}
	if (rho < 0 || rho > 1) {
		pagmo_throw(value_error,"rho must be in [0,1]");
	}
}

/// Clone method.
base_ptr aco::clone() const
{
	return base_ptr(new aco(*this));
}

/// Evolve implementation.
/**
 * Run the Firefly algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void aco::evolve(population &pop) const
{
	// Let's store some useful variables.
	problem::base_aco &prob = const_cast<problem::base_aco &>(dynamic_cast<const problem::base_aco &>(pop.problem()));
	const problem::base::size_type prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP =  pop.size();

	double max_size = 0;
	for(problem::base::size_type i = 0; i < prob_i_dimension; ++i) {
		if(max_size < (ub[i] - lb[i])) {
			max_size = ub[i] - lb[i];
		}
	}
	const std::vector<decision_vector>::size_type nComponents = boost::numeric_cast<std::vector<decision_vector>::size_type>(max_size); 

	//We perform some checks to determine wether the problem/population are suitable for ACO
	if (prob_i_dimension == 0 ) {
		pagmo_throw(value_error,"There is no integer part in the problem decision vector for Firefly to optimise");
	}

	if ( prob.get_f_dimension() != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and ACO is not suitable to solve it");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and ACO is not suitable to solve it");
	}

	if (NP < 2) {
		pagmo_throw(value_error,"for ACO at least 2 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_iter == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	decision_vector dummy(prob_i_dimension,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy);	//set of ant partial solutions
	std::vector<fitness_vector> fit(NP);		//set of ant solutions fitness

	fitness_vector tempA(prob_i_dimension,0);	//used for initialisation purpouses
	std::vector<fitness_vector> tempB(nComponents,tempA); //used for initialisation purpouses
	std::vector<std::vector<fitness_vector> > tempC(nComponents,tempB); //used for initialisation purpouses
	std::vector<std::vector<std::vector<fitness_vector> > > T(prob_i_dimension, tempC); //pheromone trail matrix 
	std::vector<std::vector<std::vector<fitness_vector> > > eta(prob_i_dimension, tempC); //pheromone trail matrix 

	// Copy the solutions and their fitness
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;
	}

	// Get heuristic information
	prob.get_heuristic_information_matrix(eta);	


	//Create pheromone paths using actual solutions
	for ( population::size_type i = 0; i<NP; i++ ) {
		deposit_pheromone(T,X[i],fit[i], m_rho);
	}


	// Main ACO loop
	for (int t = 0; t < m_iter; ++t) {

		//Select first node
		std::vector<int> selection(NP,0);
		std::vector<fitness_vector> Ttemp(nComponents*nComponents,tempA);
		for(std::vector<fitness_vector>::size_type i = 0; i < nComponents; ++i) {
			for(std::vector<fitness_vector>::size_type j = 0; j < nComponents; ++j) {
				Ttemp[i*j] = T[0][i][j];
			}
		}

		selection_probability(Ttemp, selection, pop.problem(),NP);

		for(population::size_type n=0; n < NP; ++n) {
			X[n][0] = selection[n] / nComponents + lb[0];
		}
	
		//go ahead with all the other components
		for(problem::base::size_type k = 1; k < prob_i_dimension; ++k) {
			for(population::size_type n = 0; n < NP; ++n) {
				std::vector<int> sel(1,0);
				selection_probability(T[k][X[n][k-1]], sel, pop.problem(),NP);
				X[n][k] = sel[0] + lb[k];
			}
		}
		for(population::size_type n=0; n < NP; ++n) {
			pop.set_x(n,X[n]);
			prob.objfun(fit[n], X[n]);
			deposit_pheromone(T,X[n],fit[n], m_rho);
		}
	} // end of main ACO loop

}

void aco::deposit_pheromone(std::vector<std::vector<std::vector<fitness_vector> > > &T, decision_vector &X, fitness_vector fit, double rho) {
	//evaporation
	for(std::vector<std::vector<std::vector<fitness_vector> > >::size_type k = 0; k < T.size(); ++k) {
		for(std::vector<std::vector<fitness_vector> >::size_type i=0; i < T[0].size(); ++i) {
			for(std::vector<fitness_vector>::size_type  j = 0; j < T[0][0].size(); ++j) {
				T[k][i][j][0] = (1-rho) * T[k][i][j][0];
			}
		}
	}

	//Deposit pheromone according to current solutions
	for (decision_vector::size_type i = 0; i < X.size(); ++i) {
		T[i][boost::numeric_cast<int>(X[i])][boost::numeric_cast<int>(X[i+1])][0] += rho*fit[0]; 
	}
}

//return a random integers according to the fitness probability vector in the selection vector
void aco::selection_probability(std::vector<fitness_vector> probability, std::vector<int> &selection, const pagmo::problem::base &prob, population::size_type NP) {
		std::vector<fitness_vector>::size_type pSize = probability.size();
		std::vector<double> selectionfitness(pSize), cumsum(pSize), cumsumTemp(pSize);
		double r = 0;

		fitness_vector worstfit=probability[0];
		for (std::vector<fitness_vector>::size_type i = 0; i < pSize; ++i) {
			if (prob.compare_fitness(worstfit,probability[0])) worstfit=probability[0];
		}

		for (std::vector<fitness_vector>::size_type i = 0; i < pSize; ++i) {
			selectionfitness[i] = fabs(worstfit[0] - probability[i][0]) + 1.;
		}

		// We build and normalise the cumulative sum
		cumsumTemp[0] = selectionfitness[0];
		for (std::vector<fitness_vector>::size_type i = 1; i< pSize; i++) {
			cumsumTemp[i] = cumsumTemp[i - 1] + selectionfitness[i];
		}
		for (pagmo::population::size_type i = 0; i < NP; i++) {
			cumsum[i] = cumsumTemp[i]/cumsumTemp[NP-1];
		}

		for (std::vector<double>::size_type i = 0; i < selection.size(); i++) {
			r = rng_generator::get<rng_double>()();
			for (pagmo::population::size_type j = 0; j < NP; j++) {
				if (cumsum[j] > r) {
					selection[i] = j;
					break;
				}
			}
		}
}

/// Algorithm name
std::string aco::get_name() const
{
	return "Ant Colony Optimization";
}


/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string aco::human_readable_extra() const
{
	std::ostringstream s;
	s << "iter:" << m_iter << ' ';
	s << "rho:" << m_rho << ' ';
	return s.str();
}

}} //namespaces
