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

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../problem/cstrs_self_adaptive.h"
#include "../types.h"
#include "base.h"
#include "cstrs_self_adaptive.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs an self adaptive algorithm.
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method
 * @param[in] gen number of generations.
 * @param[in] ftol stopping criteria on the f tolerance.
 * @param[in] xtol stopping criteria on the x tolerance.
 * @throws value_error if gen is negative or zero
 */
cstrs_self_adaptive::cstrs_self_adaptive(const base &original_algo, int gen, double ftol, double xtol):base(),
	m_original_algo(original_algo.clone()),m_gen(gen),m_ftol(ftol),m_xtol(xtol)
{
	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

/// Copy constructor.
cstrs_self_adaptive::cstrs_self_adaptive(const cstrs_self_adaptive &algo):base(algo),m_original_algo(algo.m_original_algo->clone()),m_gen(algo.m_gen),
	m_ftol(algo.m_ftol),m_xtol(algo.m_xtol) {}

/// Clone method.
base_ptr cstrs_self_adaptive::clone() const
{
    return base_ptr(new cstrs_self_adaptive(*this));
}

/// Evolve implementation.
/**
 * Run the Self-Adaptive algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void cstrs_self_adaptive::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type pop_size = pop.size();
	const problem::base::size_type prob_dimension = prob.get_dimension();

	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = prob.get_c_dimension();

	//We perform some checks to determine wether the problem/population are suitable for Self-Adaptive
	if(prob_c_dimension < 1) {
		pagmo_throw(value_error,"The problem is not constrained and Self-Adaptive is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (pop_size == 0) {
		return;
	}

	// Create the new problem;
    problem::cstrs_self_adaptive prob_new(prob,pop);

	// Main Self-Adaptive loop
	for(int k=0; k<m_gen; k++) {
		//std::cout << "current generation: " << k << std::endl;

		// at the first iteration the problem is not changed, 
		// for k>0 the problem is gonna change, the cache need to be reset and the population cleared
		if(k>0){ 
			prob_new.reset_caches();
			prob_new.update_penalty_coeff(pop);
		}

		// if the problem has changed it needs to be reassigned to the population
		population pop_new(prob_new,0);

		for(population::size_type i=0; i<pop_size; i++) {
			// Evaluate according to the new fitness;
			pop_new.push_back(pop.get_individual(i).cur_x);
		}

		// this constraints handling technique is initially intended to be
		// used as a fitness evaluator for ES, I am not convinced this is the easiest
		// way to implement it. But for DE, PSO in example, there is no fitness evaluator?
		m_original_algo->evolve(pop_new);

		// Reinsert best individual from the previous generation if not already present
		// Uses the constrained problem. If we don't do that, we loose the best individual...
		// Should we impose a certain percentage of the best ones?
		bool exists = false;
		for(population::size_type i=0; i<pop_size; i++) {
			if(pop_new.get_individual(i).cur_x == pop.champion().x) {
				exists=true;break;
			}
		}
		if(!exists){
			int worst=0;
			for (pagmo::population::size_type i = 1; i<pop_new.size();i++) {
				if ( prob.compare_x(pop_new.get_individual(worst).cur_x,pop_new.get_individual(i).cur_x) ) worst=i;
			}

			decision_vector dummy = pop.champion().x;
			std::transform(dummy.begin(), dummy.end(), pop.get_individual(worst).cur_x.begin(), dummy.begin(),std::minus<double>());
			//updates x and v (cache avoids to recompute the objective function)
			pop_new.set_x(worst,pop.champion().x);
			pop_new.set_v(worst,dummy);
		}

		// update the population pop
		pop.clear();
		for(pagmo::population::size_type i=0; i<pop_new.size(); i++) {
			pop.push_back(pop_new.get_individual(i).cur_x);
		}

		// Check the exit conditions (every 40 generations, just as DE)
		if(k % 40 == 0) {
			decision_vector tmp(prob_dimension);

			double dx = 0;
			for(decision_vector::size_type i=0; i<prob_dimension; i++) {
				tmp[i] = pop.get_individual(pop.get_worst_idx()).best_x[i] - pop.get_individual(pop.get_best_idx()).best_x[i];
				dx += std::fabs(tmp[i]);
			}

			if(dx < m_xtol ) {
				if (m_screen_output) {
					std::cout << "Exit condition -- xtol < " << m_xtol << std::endl;
				}
				break;
			}

			double mah = std::fabs(pop.get_individual(pop.get_worst_idx()).best_f[0] - pop.get_individual(pop.get_best_idx()).best_f[0]);

			if(mah < m_ftol) {
				if(m_screen_output) {
					std::cout << "Exit condition -- ftol < " << m_ftol << std::endl;
				}
				break;
			}

			// outputs current values
			if(m_screen_output) {
				std::cout << "Generation " << k << " ***" << std::endl;
				std::cout << "    Best global fitness: " << pop.champion().f << std::endl;
				std::cout << "    xtol: " << dx << ", ftol: " << mah << std::endl;
			}
		}

	}
}

/// Algorithm name
std::string cstrs_self_adaptive::get_name() const
{
	return m_original_algo->get_name() + "[Self-Adaptive]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr cstrs_self_adaptive::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void cstrs_self_adaptive::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string cstrs_self_adaptive::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithm: " << m_original_algo->get_name() << ' ';
	s << "\n\tConstraints handled with Self-Adaptive algorithm";
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cstrs_self_adaptive)
