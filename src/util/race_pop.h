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

#ifndef PAGMO_UTIL_RACE_POP_H
#define PAGMO_UTIL_RACE_POP_H

#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../problem/base.h"
#include "racing.h"

namespace pagmo{ namespace util {

/// racing namespace.
/**
 * Utilities for the racing mechanism.
*/
namespace racing{

/// Racing mechanism for a population
/**
 * This class implements the racing routines that can be invoked to race
 * all or some of the individuals in a population.
 *
 * It contains a caching mechanism which allows the reuse of previously
 * evaluated fitness and constraint vectors, in case an identical individual is
 * raced more than once, for example each time with different sets of other
 * individuals.  The caching mechanism ensures that all the data points that
 * are compared during the race correspond to the same seed.
 *
 * Currently the racing is implemented based on F-Race, which invokes Friedman
 * test iteratively during each race.
 *
 */
class __PAGMO_VISIBLE race_pop
{
public:

	race_pop(const population &, unsigned int seed = 0);
	race_pop(unsigned int seed = 0);

	/// Method to stop the race
	enum termination_condition { 
		 MAX_BUDGET, ///< Fixed number of function evaluations
		 MAX_DATA_COUNT ///< Fixed number of iterations
		};

	// Main method containing all the juice
	std::pair<std::vector<population::size_type>, unsigned int> run(
		const population::size_type n_final,
		const unsigned int min_trials,
		const unsigned int max_count,
		double delta,
		const std::vector<population::size_type> &,
		termination_condition term_cond,
		const bool race_best,
		const bool screen_output
	);
	
	population::size_type size() const;
	void reset_cache();
	void register_population(const population &);
	void inherit_memory(const race_pop&);
	std::vector<fitness_vector> get_mean_fitness(const std::vector<population::size_type> &active_set = std::vector<population::size_type>()) const;
	void set_seed(unsigned int);

private:
	// Helper methods to validate input data
	void _validate_active_set(const std::vector<population::size_type>& active_set, unsigned int pop_size) const;
	void _validate_problem_stochastic(const problem::base& prob) const;
	void _validate_racing_params(const population& pop, const population::size_type n_final, double delta) const;
	void _validate_budget(const unsigned int min_trials, const unsigned int max_f_evals, const std::vector<population::size_type>& in_race) const;

	unsigned int prepare_population_friedman(const std::vector<population::size_type> &in_race, unsigned int count_iter);
	unsigned int prepare_population_wilcoxon(const std::vector<population::size_type> &in_race, unsigned int count_iter);

	unsigned int compute_required_fevals(const std::vector<population::size_type>& in_race, unsigned int num_iter) const;

	// Atoms of the cache
	struct eval_data
	{
		eval_data(const fitness_vector& _f = fitness_vector(), const constraint_vector& _c = constraint_vector()): f(_f), c(_c) { }
		fitness_vector f;
		constraint_vector c;
	};

	std::vector<population::size_type> construct_output_list(
			const std::vector<racer_type>& racers,
			const std::vector<population::size_type>& decided,
			const std::vector<population::size_type>& in_race,
			const std::vector<population::size_type>& discarded,
			const population::size_type n_final,
			const bool race_best);

	// Caching routines
	void cache_insert_data(unsigned int, const fitness_vector &, const constraint_vector &);
	void cache_delete_entry(unsigned int);
	bool cache_data_exist(unsigned int, unsigned int) const;
	const eval_data &cache_get_entry(unsigned int, unsigned int) const;
	void cache_register_signatures(const population&); 
	void print_cache_stats(const std::vector<population::size_type> &) const;

	// Seeding control
	void generate_seeds(unsigned int);
	unsigned int get_current_seed(unsigned int);
	unsigned int m_race_seed;

	// Data members
	racing_population m_pop;
	racing_population m_pop_wilcoxon;
	bool m_pop_registered;
	std::vector<unsigned int> m_seeds;
	rng_uint32 m_seeder;
	bool m_use_caching;
	std::vector<std::vector<eval_data> > m_cache_data;
	std::vector<eval_data> m_cache_averaged_data;
	std::vector<decision_vector> m_cache_signatures;
};

}}}

#endif
