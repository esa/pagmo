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
#include <sstream>
#include <string>

#include "algorithm/base.h"
#include "problem/base.h"
#include "exceptions.h"
#include "island.h"
#include "island_storage.h"
#include "types.h"

namespace pagmo
{

/// Constructor from problem::base, algorithm::base and number of individuals.
/**
 * Will store a copy of the problem and of the algorithm internally, will initialise internal population to n individuals
 * and evolution time to zero. Will fail if n is negative.
 */
island::island(const problem::base &p, const algorithm::base &a, int n):island_storage(p,a,n),m_archi(0),m_evo_time(0) {}

/// Copy constructor.
/**
 * Will perform a deep copy of all the elements of island isl, which will be synchronised before any operation takes place.
 */
island::island(const island &isl):island_storage()
{
	operator=(isl);
}

/// Assignment operator.
/**
 * Performs a deep copy of all the elements of isl into this island. Both island will be synchronised before assignment.
 */
island &island::operator=(const island &isl)
{
	if (this != &isl) {
		// Make sure both islands are in a known state.
		join();
		isl.join();
		// Call the operator from island_storage.
		island_storage::operator=(isl);
		// Archipelago pointer.
		m_archi = isl.m_archi;
		// Copy over also evolution time.
		m_evo_time = isl.m_evo_time;
	}
	return *this;
}

/// Destructor.
/**
 * Will call island::join() before returning.
 */
island::~island()
{
	join();
}

/// Explicitly synchronise the island.
/**
 * After this method returns, any pending evolution has been completed.
 */
void island::join() const
{
	lock_type lock(m_evo_mutex);
}

/// Return terse human readable representation of the island.
/**
 * Will return a formatted string containing the human readable forms of the problem and of the algorithm.
 */
std::string island::human_readable_terse() const
{
	join();
	std::ostringstream oss;
	oss << prob() << '\n';
	oss << algo() << '\n';
	return oss.str();
}

/// Return complete human readable representation of the island.
/**
 * Will return a formatted string containing the output of human_readable_terse(),
 * whether or not the island belongs to an archipelago, the total evolution time and the list of individuals.
 */
std::string island::human_readable() const
{
	join();
	std::ostringstream oss;
	oss << human_readable_terse();
	oss << "Belongs to archipelago: " << (m_archi ? "true" : "false") << '\n' << '\n';
	if (pop().size()) {
		oss << "List of individuals:\n";
		for (size_type i = 0; i < pop().size(); ++i) {
			oss << '#' << i << ":\n";
			oss << "\tDecision vector:\t" << pop()[i].get<0>() << '\n';
			oss << "\tVelocity vector:\t" << pop()[i].get<1>() << '\n';
			oss << "\tFitness vector:\t\t" << pop()[i].get<2>() << '\n';
			oss << "\tBest fitness vector:\t" << pop()[i].get<3>() << '\n';
		}
	} else {
		oss << "No individuals.\n";
	}
	return oss.str();
}

/// Overload stream operator for pagmo::island.
/**
 * Equivalent to printing island::human_readable() to stream.
 */
std::ostream &operator<<(std::ostream &s, const island &isl)
{
	s << isl.human_readable();
	return s;
}

}
