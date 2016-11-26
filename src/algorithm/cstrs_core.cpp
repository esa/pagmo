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
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../problem/con2uncon.h"
#include "../types.h"
#include "base.h"
#include "cstrs_core.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs a CORE constraints handling algorithm
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method. Its number of
 * generations should be set to 1.
 * @param[in] repair_algo pagmo::algorithm to use as 'repairing' optimization algorithm. It should be
 * able to deal with population of size 1.
 * @param[in] gen number of generations.
 * @param[in] repair_frequency The infeasible are repaired at each repair frequency generations.
 * @param[in] repair_ratio It the repair ratio is the ratio of repaired individuals over infeasible
 * ones (a ratio of 1 will repair all the individuals).
 * @param[in] ftol stopping criteria on the f tolerance.
 * @param[in] xtol stopping criteria on the x tolerance.
 * @throws value_error if gen is negative, if repair frequency is negative.
 */
cstrs_core::cstrs_core(const base &original_algo, const base &repair_algo,
                       int gen,
					   int repair_frequency,
					   double repair_ratio,
					   double ftol, double xtol):
    base(),m_original_algo(original_algo.clone()),
    m_repair_algo(repair_algo.clone()),
    m_gen(gen),m_repair_frequency(repair_frequency),
	m_repair_ratio(repair_ratio),m_ftol(ftol),m_xtol(xtol)
{
//	m_original_algo = original_algo.clone();
//    m_repair_algo = repair_algo.clone();

	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if(repair_frequency < 0) {
		pagmo_throw(value_error,"repair frequency must be positive");
	}
	if((repair_ratio < 0) || (repair_ratio > 1)) {
		pagmo_throw(value_error,"repair ratio must be in [0..1]");
	}
}

/// Copy constructor.
cstrs_core::cstrs_core(const cstrs_core &algo):
    base(algo),m_original_algo(algo.m_original_algo->clone()),
    m_repair_algo(algo.m_repair_algo->clone()),m_gen(algo.m_gen),
    m_repair_frequency(algo.m_repair_frequency),
    m_repair_ratio(algo.m_repair_ratio),
    m_ftol(algo.m_ftol),m_xtol(algo.m_xtol)
{}

/// Clone method.
base_ptr cstrs_core::clone() const
{
	return base_ptr(new cstrs_core(*this));
}

/// Evolve implementation.
/**
 * Run the CORE algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void cstrs_core::evolve(population &pop) const
{	
	// store useful variables
	const problem::base &prob = pop.problem();
	const population::size_type pop_size = pop.size();
	const problem::base::size_type prob_dimension = prob.get_dimension();

	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = prob.get_c_dimension();

	//We perform some checks to determine wether the problem/population are suitable for CORE
	if(prob_c_dimension < 1) {
		pagmo_throw(value_error,"The problem is not constrained and CORE is not suitable to solve it");
	}
	if(prob.get_f_dimension() != 1) {
		pagmo_throw(value_error,"The problem is multiobjective and CORE is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if(pop_size == 0) {
		return;
	}

	// generates the unconstrained problem
	problem::con2uncon prob_unconstrained(prob);

	// associates the population to this problem
	population pop_uncon(prob_unconstrained);

	// fill this unconstrained population
	pop_uncon.clear();
	for(population::size_type i=0; i<pop_size; i++) {
		pop_uncon.push_back(pop.get_individual(i).cur_x);
	}

	// vector containing the infeasibles positions
	std::vector<population::size_type> pop_infeasibles;

	// Main CORE loop
	for(int k=0; k<m_gen; k++) {

		if(k%m_repair_frequency == 0) {
			pop_infeasibles.clear();

			// get the infeasible individuals
			for(population::size_type i=0; i<pop_size; i++) {
				if(!prob.feasibility_c(pop.get_individual(i).cur_c)) {
					pop_infeasibles.push_back(i);
				}
			}

			// random shuffle of infeasibles?
			population::size_type number_of_repair = (population::size_type)(m_repair_ratio * pop_infeasibles.size());

			// repair the infeasible individuals
			for(population::size_type i=0; i<number_of_repair; i++) {
				const population::size_type &current_individual_idx = pop_infeasibles.at(i);

                pop.repair(current_individual_idx, m_repair_algo);
			}

			// the population is repaired, it can be now used in the new unconstrained population
			// only the repaired individuals are put back in the population
			for(population::size_type i=0; i<number_of_repair; i++) {
				population::size_type current_individual_idx = pop_infeasibles.at(i);
				pop_uncon.set_x(current_individual_idx, pop.get_individual(current_individual_idx).cur_x);
			}
		}

		m_original_algo->evolve(pop_uncon);

		// push back the population in the main problem
		pop.clear();
		for(population::size_type i=0; i<pop_size; i++) {
			pop.push_back(pop_uncon.get_individual(i).cur_x);
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
std::string cstrs_core::get_name() const
{
	return m_original_algo->get_name() + "[CORE]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr cstrs_core::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void cstrs_core::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Get a copy of the internal local repair algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local repair algorithm.
 */
base_ptr cstrs_core::get_repair_algorithm() const
{
    return m_repair_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local repair algorithm.
 *
 * @param[in] algo algorithm to be set as local repair algorithm.
 */
void cstrs_core::set_repair_algorithm(const base &algo)
{
    m_repair_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string cstrs_core::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithms: " << m_original_algo->get_name() << " ";
	s << "\n\tConstraints handled with CORE algorithm";
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cstrs_core)
