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
#include "../types.h"
#include "base.h"
#include "ms.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] algorithm pagmo::algorithm for the multistarts
 * @param[in] starts number of multistarts
 * @throws value_error if starts is negative
 */
ms::ms(const base &algorithm, int starts):base(),m_starts(starts)
{
	m_algorithm = algorithm.clone();
	if (starts < 0) {
		pagmo_throw(value_error,"number of multistarts needs to be larger than zero");
	}
}

/// Copy constructor (deep copy).
ms::ms(const ms &other):base(other),m_algorithm(other.m_algorithm->clone()),m_starts(other.m_starts) {}

/// Clone method.
base_ptr ms::clone() const
{
	return base_ptr(new ms(*this));
}

/// Evolve implementation.
/**
 * Run the Multi-start algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void ms::evolve(population &pop) const
{
	// Let's store some useful variables.
	const population::size_type NP = pop.size();

	// Get out if there is nothing to do.
	if (m_starts == 0 || NP == 0) {
		return;
	}

	// Local population used in the algorithm iterations.
	population working_pop(pop);

	//ms main loop
	for (int i=0; i< m_starts; ++i)
	{
		working_pop.reinit();
		m_algorithm->evolve(working_pop);
		if (working_pop.problem().compare_fc(working_pop.get_individual(working_pop.get_best_idx()).cur_f,working_pop.get_individual(working_pop.get_best_idx()).cur_c,
			pop.get_individual(pop.get_worst_idx()).cur_f,pop.get_individual(pop.get_worst_idx()).cur_c
		) )
		{
			//update best population replacing its worst individual with the good one just produced.
			pop.set_x(pop.get_worst_idx(),working_pop.get_individual(working_pop.get_best_idx()).cur_x);
			pop.set_v(pop.get_worst_idx(),working_pop.get_individual(working_pop.get_best_idx()).cur_v);
		}
		if (m_screen_output)
		{
			std::cout << i << ". " << "\tCurrent iteration best: " << working_pop.get_individual(working_pop.get_best_idx()).cur_f << "\tOverall champion: " << pop.champion().f << std::endl;
		}
	}
}


/// Algorithm name
std::string ms::get_name() const
{
	return "Multi-start";
}

/// Get a copy of the internal algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal algorithm.
 */
base_ptr ms::get_algorithm() const
{
	return m_algorithm->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal algorithm.
 * 
 * @param[in] algo algorithm to be set for multistart.
 */
void ms::set_algorithm(const base &algo)
{
	m_algorithm = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string ms::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithm: " << m_algorithm->get_name() << ' ';
	s << "iter:" << m_starts << ' ';
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::ms)
