/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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
 * This class contains ...
 */
class __PAGMO_VISIBLE race_pop
{
public:

	race_pop(const population&, unsigned int seed = 0);

	// Main method containing all the juice
	std::pair<std::vector<population::size_type>, unsigned int>  run(
		const population::size_type n_final,
		const unsigned int min_trials,
		const unsigned int max_count,
		double delta,
		const std::vector<population::size_type> &,
		const bool race_best,
		const bool screen_output
	);

	void reset_cache();

private:
	// Helper methods to validate input data
	void _validate_active_set(const std::vector<population::size_type>& active_set, unsigned int pop_size);
	void _validate_problem_stochastic(const problem::base& prob);
	void _validate_racing_params(const population& pop, const population::size_type n_final, const unsigned int min_trials, const unsigned int max_f_evals, double delta, unsigned int active_set_size);

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

	// Seeding control
	void generate_seeds(unsigned int);
	unsigned int get_current_seed(unsigned int);
	
	// Data members
	racing_population m_pop;
	std::vector<unsigned int> m_seeds;
	rng_uint32 m_seeder;
	std::vector<std::vector<eval_data> > m_cache_data;
};

}}}

#endif
