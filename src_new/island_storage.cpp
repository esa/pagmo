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

// 04/01/2009: Initial version by Francesco Biscani.

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>

#include "algorithm/base.h"
#include "problem/base.h"
#include "exceptions.h"
#include "island_storage.h"
#include "rng.h"
#include "types.h"

namespace pagmo
{
/// Default constructor.
/**
 * Will default-initialise all data members.
 */
island_storage::island_storage() {}

// Function object to compare individuals according to their current fitness.
struct cur_f_comp {
	cur_f_comp(const problem::base &p):m_p(p) {}
	bool operator()(const island_storage::individual_type &i1, const island_storage::individual_type &i2) const
	{
		return m_p.compare_f(i1.get<2>(),i2.get<2>());
	}
	const problem::base &m_p;
};

/// Constructor from problem::base, algorithm::base and number of individuals.
/**
 * Will store a copy of the problem and of the algorithm internally, and it will initialise the internal population to n randomly-generated individuals.
 * Will fail if n is negative.
 */
island_storage::island_storage(const problem::base &p, const algorithm::base &a, int n):m_prob(p.clone()),m_algo(a.clone()),m_drng(rng_generator::get<rng_double>())
{
	if (n < 0) {
		pagmo_throw(value_error,"number of individuals cannot be negative");
	}
	// Store sizes temporarily.
	const size_type size = boost::numeric_cast<size_type>(n);
	const fitness_vector::size_type f_size = m_prob->get_f_dimension();
	const decision_vector::size_type p_size = m_prob->get_dimension(), i_size = m_prob->get_i_dimension();
	pagmo_assert(p_size >= i_size);
	for (size_type i = 0; i < size; ++i) {
		// Push back an empty individual.
		m_pop.push_back(individual_type());
		// Resize decision vectors and fitness vectors.
		m_pop.back().get<0>().resize(p_size);
		m_pop.back().get<1>().resize(p_size);
		m_pop.back().get<2>().resize(f_size);
		m_pop.back().get<3>().resize(f_size);
		// Initialise randomly the continuous part of the decision vector.
		for (decision_vector::size_type i = 0; i < p_size - i_size; ++i) {
			m_pop.back().get<0>()[i] = m_prob->get_lb()[i] + m_drng() * (m_prob->get_ub()[i] - m_prob->get_lb()[i]);
		}
		// Initialise randomly the integer part of the decision vector.
		for (decision_vector::size_type i = p_size - i_size; i < p_size; ++i) {
			m_pop.back().get<0>()[i] = double_to_int::nearbyint(m_prob->get_lb()[i] + m_drng() * (m_prob->get_ub()[i] - m_prob->get_lb()[i]));
		}
		// Initialise randomly the velocity vector.
		for (decision_vector::size_type i = 0; i < p_size; ++i) {
			m_pop.back().get<1>()[i] = (m_drng() - .5) * (m_prob->get_ub()[i] - m_prob->get_lb()[i]);
		}
		// Compute the current fitness.
		m_prob->objfun(m_pop.back().get<2>(),m_pop.back().get<0>());
		// Best fitness is current fitness.
		m_pop.back().get<3>() = m_pop.back().get<2>();
	}
	// Calculate the champion.
	population_type::iterator it = std::min_element(m_pop.begin(),m_pop.end(),cur_f_comp(*m_prob));
	if (it != m_pop.end()) {
		m_champion.get<0>() = it->get<0>();
		m_champion.get<1>() = it->get<2>();
	}
}

/// Assignment operator.
/**
 * Performs a deep copy of all the elements of isl into this.
 */
island_storage &island_storage::operator=(const island_storage &isl)
{
	if (this != &isl) {
		// Perform the copies.
		m_prob = isl.m_prob->clone();
		m_algo = isl.m_algo->clone();
		m_pop = isl.m_pop;
		m_champion = isl.m_champion;
		m_drng = isl.m_drng;
	}
	return *this;
}

/// Get constant reference to internal problem::base class.
const problem::base &island_storage::prob() const
{
	return *m_prob;
}

/// Get constant reference to internal algorithm::base class.
const algorithm::base &island_storage::algo() const
{
	return *m_algo;
}

/// Get constant reference to internal population.
const island_storage::population_type &island_storage::pop() const
{
	return m_pop;
}

/// Get constant reference to internal champion member.
const island_storage::champion_type &island_storage::champion() const
{
	return m_champion;
}

}
