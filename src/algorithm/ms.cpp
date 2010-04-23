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
ms::ms(const algorithm::base &algorithm, int starts):base(),m_starts(starts)
{
	m_algorithm = algorithm.clone();
	if (starts < 0) {
		pagmo_throw(value_error,"number of multistarts needs to be larger than zero");
	}
}

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

	// Init the best fitness and constraint vector
	fitness_vector best_f = pop.get_individual(pop.get_best_idx()).cur_f;
	constraint_vector best_c = pop.get_individual(pop.get_best_idx()).cur_c;
	population best_pop(pop);

	//ms main loop
	for (int i=0; i< m_starts; ++i)
	{
		pop.reinit();
		m_algorithm->evolve(pop);
		if (m_screen_out)
		{
			std::cout << i << ". " << "\tBest: " << pop.get_individual(pop.get_best_idx()).cur_f << "\tChampion: " << pop.champion().f << std::endl;
		}

		if (pop.problem().compare_fc(pop.get_individual(pop.get_best_idx()).cur_f,pop.get_individual(pop.get_best_idx()).cur_c,best_f,best_c) )
		{
			best_f = pop.get_individual(pop.get_best_idx()).cur_f;
			best_c = pop.get_individual(pop.get_best_idx()).cur_c;

			//update best population
			for (population::size_type j=0; j<pop.size();++j)
			{
				best_pop.set_x(j,pop.get_individual(j).cur_x);
				best_pop.set_v(j,pop.get_individual(j).cur_v);
			}

		}
	}
	//on exit set the population to the best one of the multistarts
	for (population::size_type j=0; j<pop.size();++j)
	{
		pop.set_x(j,best_pop.get_individual(j).cur_x);
		pop.set_v(j,best_pop.get_individual(j).cur_v);
	}
}

/// Activate screen output
/**
 * Activate screen output. Everytime a new champion is found the following information is printed
 * on the screen: Number of iterations, best fitness
 *
 * @param[in] p true or false
 */
void ms::screen_output(const bool p) {m_screen_out = p;}


/// Algorithm name
std::string ms::get_name() const
{
	return "Multi-start";
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string ms::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithm:" << m_algorithm->get_name() << ' ';
	s << "starts:" << m_starts << ' ';
	return s.str();
}

}} //namespaces
