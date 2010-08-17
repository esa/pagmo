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
#include "../problem/base_aco.h"
#include "../types.h"




namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] iter number of iterations.
 * @param[in] rho evaporation rate
 * @throws value_error if number of iterations is negative or rho isn't in the [0,1] range.
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
 * Run the ACO algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void aco::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base_aco &prob = dynamic_cast<const problem::base_aco &>(pop.problem());
	const problem::base::size_type prob_i_dimension = prob.get_i_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP =  pop.size();

	double max_size = 0;
	for(problem::base::size_type i = 0; i < prob_i_dimension; ++i) {
		if(max_size < (ub[i] - lb[i])) {
			max_size = ub[i] - lb[i];
		}
	}
	const std::vector<decision_vector>::size_type nComponents = boost::numeric_cast<std::vector<decision_vector>::size_type>(max_size) + 1; 

	//We perform some checks to determine wether the problem/population are suitable for ACO
	if (prob_i_dimension == 0 ) {
		pagmo_throw(value_error,"There is no integer part in the problem decision vector for ACO to optimise");
	}

	if ( prob.get_f_dimension() != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and ACO is not suitable to solve it");
	}

	if (NP < 2) {
		pagmo_throw(value_error,"for ACO at least 2 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_iter == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	decision_vector dummy(prob_i_dimension,0);	//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy);	//set of ant partial solutions
	std::vector<fitness_vector> fit(NP);		//set of ant solutions fitness

	fitness_vector tempA(prob.get_f_dimension(),0);	//used for initialisation purpouses
	std::vector<fitness_vector> tempB(nComponents,tempA); //used for initialisation purpouses
	std::vector<std::vector<fitness_vector> > tempC(nComponents,tempB); //used for initialisation purpouses
	std::vector<std::vector<std::vector<fitness_vector> > > T(prob_i_dimension, tempC); //pheromone trail matrix 
	std::vector<std::vector<std::vector<fitness_vector> > > eta = prob.get_heuristic_information_matrix(); //heuristic information matrix 

	// Copy the solutions and their fitness
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;
	}

	//Create pheromone paths using actual solutions
	for ( population::size_type i = 0; i<NP; i++ ) {
		deposit_pheromone(T,X[i],fit[i], m_rho);
	}


	std::vector<int> selection(NP,0); //next node selection for each individual
	std::vector<fitness_vector> Ttemp(nComponents,tempA);
	std::vector<fitness_vector> etaTemp(nComponents,tempA);
	std::vector<bool> fComponentsTemp(nComponents,true);
	std::vector<bool> fComponents(nComponents);

	// Main ACO loop
	for (int t = 0; t < m_iter; ++t) {

		//Select first node
		// Since the first node doesn't have a predecessor we create new T and eta vector for each first vector i
		// summing over all the j each T[i][j] 
		for(std::vector<fitness_vector>::size_type i = 0; i < nComponents; ++i) {
			for(std::vector<fitness_vector>::size_type j = 0; j < nComponents; ++j) {
				Ttemp[i][0] += T[0][i][j][0];
				etaTemp[i][0] += eta[0][i][j][0];
			}
		}

		selection_probability(Ttemp, fComponentsTemp, etaTemp, selection, pop.problem());

		for(population::size_type n=0; n < NP; ++n) {
			X[n][0] = selection[n] + lb[0];
		}
	
		//go ahead with all the other components
		for(problem::base::size_type k = 1; k < prob_i_dimension; ++k) {
			for(population::size_type n = 0; n < NP; ++n) {
				feasible_components(fComponents, prob, X[n], k, lb[k], ub[k]);
				std::vector<int> sel(1,0);
				selection_probability(T[k][X[n][k-1]], fComponents, eta[k][X[n][k-1]], sel, pop.problem());
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

/*
 * given a partial decision vector of size X.size() that has just the first xSize components filled (it means that the rest of the vector is not relevant)
 * write a fComonents boolean vector where fComponents[i] is true if X[xSize] =  "ith possible value for the component xSize" make
 * the vector X feasible (according to the first xSize+1 components)
 */
void aco::feasible_components(std::vector<bool> &fComponents, const pagmo::problem::base_aco &prob, decision_vector &X, problem::base::size_type xSize, double lb, double ub) {
	decision_vector tmpX(xSize+1,0);
	for(pagmo::problem::base::size_type i = 0; i < xSize; ++i) {
		tmpX[i] = X[i];
	}
	int i = 0;
	for (int n = boost::numeric_cast<int>(lb); i+n <= ub; ++i) {
		tmpX[xSize] = i+n;
		if (prob.check_partial_feasibility(tmpX)) {
			fComponents[i] = true;
		}
		else {
			fComponents[i] = false;
		}
	}
}

/*
 * Deposit pherormone on the trail. Pheromone is represented as the fitness of a solution. Each individual deposit an amount of pheromone
 * on its path (its solution) equal to the fitness of its solution
 */
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
	T[0][boost::numeric_cast<int>(X[0])][boost::numeric_cast<int>(X[1])][0] += rho*fit[0]; 
	for (decision_vector::size_type i = 1; i < X.size(); ++i) {
		T[i-1][boost::numeric_cast<int>(X[i-1])][boost::numeric_cast<int>(X[i])][0] += rho*fit[0]; 
	}
}

//return a random integers according to the fitness probability vector in the selection vector
void aco::selection_probability(std::vector<fitness_vector> &probability, std::vector<bool> &fComponents, std::vector<fitness_vector> &eta, std::vector<int> &selection, const pagmo::problem::base &prob) {
		std::vector<fitness_vector>::size_type pSize = probability.size();
		std::vector<double> selectionfitness(pSize), cumsum(pSize), cumsumTemp(pSize);
		double r = 0;

		fitness_vector worstfit(1, probability[0][0]);
		fitness_vector tmpFit(1, probability[0][0]);
		for (std::vector<fitness_vector>::size_type i = 0; i < pSize; ++i) {
			tmpFit[0] = probability[i][0] * eta[i][0];
			if (prob.compare_fitness(worstfit,tmpFit)) worstfit[0]=tmpFit[0];
		}

		for (std::vector<fitness_vector>::size_type i = 0; i < pSize; ++i) {
			if(fComponents[i]) {
				selectionfitness[i] = fabs(worstfit[0] - eta[i][0]*probability[i][0]) + 1.;
			}
			else {
				selectionfitness[i] = 0;
			}
		}

		// We build and normalise the cumulative sum
		cumsumTemp[0] = selectionfitness[0];
		for (std::vector<fitness_vector>::size_type i = 1; i< pSize; i++) {
			cumsumTemp[i] = cumsumTemp[i - 1] + selectionfitness[i];
		}
		for (pagmo::population::size_type i = 0; i < pSize; i++) {
			cumsum[i] = cumsumTemp[i]/cumsumTemp[pSize-1];
		}

		for (std::vector<double>::size_type i = 0; i < selection.size(); i++) {
			r = rng_generator::get<rng_double>()();
			for (pagmo::population::size_type j = 0; j < pSize; j++) {
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
