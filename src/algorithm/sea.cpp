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
#include <string>
#include <vector>
#include <algorithm>

#include "../exceptions.h"
#include "../population.h"
#include "../types.h"
#include "base.h"
#include "sea.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen Number of generations to evolve.
 * @throws value_error if gen is negative, or the problem i not integer, box contrained and single-objective
 *
 */
sea::sea(int gen)
	:base(),m_gen(gen)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

/// Clone method.
base_ptr sea::clone() const
{
	return base_ptr(new sea(*this));
}

/// Evolve implementation.
/**
 * Run the EA for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void sea::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type Di = prob.get_i_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();

	//We perform some checks to determine wether the problem/population are suitable for SGA
	if ( prob.get_c_dimension() != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and EA is not suitable to solve it");
	}

	if ( prob.get_dimension() - prob.get_i_dimension() != 0 ) {
		pagmo_throw(value_error,"The problem has a continuous dimension and this (N+1)-EA Simple Evolutionary Algorithm is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}

	int new_gene;
	// Main loop
	for (int j = 0; j<m_gen; j++) {

		// Offspring is generated from the best individual
		decision_vector offspring = pop.get_individual(pop.get_best_idx()).cur_x;

		// Mutation of the best individual. Each gene is flipped with probability 1/Di to a different value
		for (pagmo::problem::base::size_type j = 0; j < Di;++j) {//for each integer variable
			if (m_drng() < 1.0/Di) {
				do {
					new_gene = boost::uniform_int<int>(lb[j],ub[j])(m_urng);
				} while(new_gene == offspring[j]);
				offspring[j] = new_gene;
			}
		}
		
		// We add the mutated individual to the population (this will also evaluate its fitness)
		pop.push_back(offspring);
		// We get rid of the worst individual (in multi-objective this is computed using
		// the crowding distance operator)
		pop.erase(pop.get_worst_idx());

	} // end of main loop
}

/// Algorithm name
std::string sea::get_name() const
{
	return "(N+1)-EA Simple Evolutionary Algorithm";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string sea::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::sea)
