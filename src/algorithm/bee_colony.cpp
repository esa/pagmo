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

#include "bee_colony.h"
#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"




namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations (2 * pop.size() function evaluations per generation).
 * @param[in] limit number of tries after which a source of food is dropped if not improved
 * @throws value_error if number of iterations or limit are negative
 */
bee_colony::bee_colony(int gen, int limit):base(),m_iter(gen), m_limit(limit) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}

	if (limit < 0) {
		pagmo_throw(value_error,"limit value must be nonnegative");
	}

}

/// Clone method.
base_ptr bee_colony::clone() const
{
	return base_ptr(new bee_colony(*this));
}

/// Evolve implementation.
/**
 * Run the ABC algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void bee_colony::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type prob_i_dimension = prob.get_i_dimension(), D = prob.get_dimension(), Dc = D - prob_i_dimension, prob_c_dimension = prob.get_c_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = (int) pop.size();

	//We perform some checks to determine whether the problem/population are suitable for ABC
	if ( Dc == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for ABC to optimise");
	}

	if ( prob.get_f_dimension() != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and ABC is not suitable to solve it");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and ABC is not suitable to solve it");
	}

	if (NP < 2) {
		pagmo_throw(value_error,"for ABC at least 2 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_iter == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	fitness_vector fnew(prob.get_f_dimension());
	decision_vector dummy(D,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy);	//set of food sources
	std::vector<fitness_vector> fit(NP);		//food sources fitness

	decision_vector temp_solution(D,0);

	std::vector<int> trial(NP,0);

	std::vector<double> probability(NP);

	population::size_type neighbour = 0;

	decision_vector::size_type param2change = 0;

	std::vector<double> selectionfitness(NP), cumsum(NP), cumsumTemp(NP);
	std::vector <population::size_type> selection(NP);


	double r = 0;

	// Copy the food sources position and their fitness
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;
	}

	// Main ABC loop
	for (int j = 0; j < m_iter; ++j) {
		//1- Send employed bees
		for (population::size_type ii = 0; ii< NP; ++ii) {
			//selects a random component (only of the continuous part) of the decision vector
			param2change = boost::uniform_int<decision_vector::size_type>(0,Dc-1)(m_urng);
			//randomly chose a solution to be used to produce a mutant solution of solution ii
			//randomly selected solution must be different from ii
			do{
				neighbour = boost::uniform_int<population::size_type>(0,NP-1)(m_urng);
			}
			while(neighbour == ii);

			//copy local solution into temp_solution (the whole decision_vector, also the integer part)
			for(population::size_type i=0; i<D; ++i) {
				temp_solution[i] = X[ii][i];
			}

			//mutate temp_solution
			temp_solution[param2change] = X[ii][param2change] + boost::uniform_real<double>(-1,1)(m_drng) * (X[ii][param2change] - X[neighbour][param2change]);

			//if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
			if (temp_solution[param2change]<lb[param2change]) {
				temp_solution[param2change] = lb[param2change];
			}
			if (temp_solution[param2change]>ub[param2change]) {
				temp_solution[param2change] = ub[param2change];
			}

			//Calling void prob.objfun(fitness_vector,decision_vector) is more efficient as no memory allocation occur
			//A call to fitness_vector prob.objfun(decision_vector) allocates memory for the return value.
			prob.objfun(fnew,temp_solution);
			//If the new solution is better than the old one replace it with the mutant one and reset its trial counter
			if(prob.compare_fitness(fnew, fit[ii])) {
				X[ii][param2change] = temp_solution[param2change];
				pop.set_x(ii,X[ii]);
				prob.objfun(fit[ii], X[ii]); //update the fitness vector
				trial[ii] = 0;
			}
			else {
				trial[ii]++; //if the solution can't be improved increase its trial counter
			}
		} //End of loop on the population members

		//2 - Send onlooker bees
		//We scale all fitness values from 0 (worst) to absolute value of the best fitness
		fitness_vector worstfit=fit[0];
		for (pagmo::population::size_type i = 1; i < NP;i++) {
			if (prob.compare_fitness(worstfit,fit[i])) worstfit=fit[i];
		}

		for (pagmo::population::size_type i = 0; i < NP; i++) {
			selectionfitness[i] = fabs(worstfit[0] - fit[i][0]) + 1.;
		}

		// We build and normalise the cumulative sum
		cumsumTemp[0] = selectionfitness[0];
		for (pagmo::population::size_type i = 1; i< NP; i++) {
			cumsumTemp[i] = cumsumTemp[i - 1] + selectionfitness[i];
		}
		for (pagmo::population::size_type i = 0; i < NP; i++) {
			cumsum[i] = cumsumTemp[i]/cumsumTemp[NP-1];
		}

		for (pagmo::population::size_type i = 0; i < NP; i++) {
			r = m_drng();
			for (pagmo::population::size_type j = 0; j < NP; j++) {
				if (cumsum[j] > r) {
					selection[i]=j;
					break;
				}
			}
		}

		for(pagmo::population::size_type t = 0; t < NP; ++t) {
			r = m_drng();
			pagmo::population::size_type ii = selection[t];
			//selects a random component (only of the continuous part) of the decision vector
			param2change = boost::uniform_int<decision_vector::size_type>(0,Dc-1)(m_urng);
			//randomly chose a solution to be used to produce a mutant solution of solution ii
			//randomly selected solution must be different from ii
			do{
				neighbour = boost::uniform_int<population::size_type>(0,NP-1)(m_urng);
			}
			while(neighbour == ii);

			//copy local solution into temp_solution (also integer part)
			for(population::size_type i=0; i<D; ++i) {
				temp_solution[i] = X[ii][i];
			}

			//mutate temp_solution
			temp_solution[param2change] = X[ii][param2change] + boost::uniform_real<double>(-1,1)(m_drng) * (X[ii][param2change] - X[neighbour][param2change]);

			/*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
			if (temp_solution[param2change]<lb[param2change]) {
				temp_solution[param2change] = lb[param2change];
			}
			if (temp_solution[param2change]>ub[param2change]) {
				temp_solution[param2change] = ub[param2change];
			}

			//Calling void prob.objfun(fitness_vector,decision_vector) is more efficient as no memory allocation occur
			//A call to fitness_vector prob.objfun(decision_vector) allocates memory for the return value.
			prob.objfun(fnew,temp_solution);
			//If the new solution is better than the old one replace it with the mutant one and reset its trial counter
			if(prob.compare_fitness(fnew, fit[ii])) {
				X[ii][param2change] = temp_solution[param2change];
				pop.set_x(ii,X[ii]);
				prob.objfun(fit[ii], X[ii]); //update the fitness vector
				trial[ii] = 0;
			}
			else {
				trial[ii]++; //if the solution can't be improved increase its trial counter
			}
		}

		//3 - Send scout bees
		int maxtrialindex = 0;
		for (population::size_type ii=1; ii<NP; ++ii)
		{
			if (trial[ii] > trial[maxtrialindex]) {
				maxtrialindex = ii;
			}
		}
		if(trial[maxtrialindex] >= m_limit)
		{
			//select a new random solution
			for(problem::base::size_type jj = 0; jj < Dc; ++jj) {
				X[maxtrialindex][jj] = boost::uniform_real<double>(lb[jj],ub[jj])(m_drng);
			}
			trial[maxtrialindex] = 0;
			pop.set_x(maxtrialindex,X[maxtrialindex]);
		}

	} // end of main ABC loop

}

/// Algorithm name
std::string bee_colony::get_name() const
{
	return "Artificial Bee Colony optimization";
}


/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string bee_colony::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_iter << ' ';
	s << "limit:" << m_limit << ' ';
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::bee_colony)
