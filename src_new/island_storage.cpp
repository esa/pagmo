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

#include <boost/numeric/conversion/cast.hpp>

#include "algorithm/base.h"
#include "problem/base.h"
#include "exceptions.h"
#include "island_storage.h"
#include "types.h"

namespace pagmo
{
/// Default constructor.
/**
 * Not intended for stand-alone use.
 */
island_storage::island_storage() {}


/// Constructor from problem::base, algorithm::base and number of individuals.
/**
 * Will store a copy of the problem and of the algorithm internally, will initialise internal population to n individuals
 * and evolution time to zero. Will fail if n is negative.
 */
island_storage::island_storage(const problem::base &p, const algorithm::base &a, int n):m_prob(p.clone()),m_algo(a.clone())
{
	if (n < 0) {
		pagmo_throw(value_error,"number of individuals cannot be negative");
	}
	const size_type size = boost::numeric_cast<size_type>(n);
	for (size_type i = 0; i < size; ++i) {
		
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

}
