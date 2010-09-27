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

#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>

#include "../population.h"
#include "../rng.h"
#include "base.h"

namespace pagmo
{
namespace algorithm {

/// Default constructor.
/**
 * Will initialise the internal random number generators using seeds provided by the thread-safe pagmo::static_rng_uint32.
 */
base::base():m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>()) {}

/// Trivial destructor.
/**
 * No side effects.
 */
base::~base() {}

/// Get algorithm's name.
/**
 * Default implementation will return the algorithm's mangled C++ name.
 *
 * @return name of the algorithm.
 */
std::string base::get_name() const
{
	return typeid(*this).name();
}

/// Return human readable representation of the algorithm.
/**
 * Will return a formatted string containing the algorithm name from get_name().
 * The output of human_readable_extra() will be appended at the end of the string.
 *
 * @return string containing human readable representation of the algorithm.
 */
std::string base::human_readable() const
{
	std::ostringstream s;
	s << get_name() << " - ";
	s << human_readable_extra();
	return s.str();
}

/// Extra information in human readable format.
/**
 * Default implementation returns an empty string.
 *
 * @return string containing algorithm specific human readable representation of the algorithm.
 */
std::string base::human_readable_extra() const
{
	return std::string();
}

/// Algorithm's thread safety property.
/**
 * Return true if the algorithm is thread-safe.
 * Default implementation returns true.
 *
 * @return true if the algorithm is thread-safe, false otherwise.
 */
bool base::is_thread_safe() const
{
	return true;
}

/// Overload stream operator for algorithm::base.
/**
 * Equivalent to printing base::human_readable() to stream.
 *
 * @param[in] s stream to which the algorithm will be sent.
 * @param[in] alg algorithm::base to be sent to stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const base &alg)
{
	s << alg.human_readable();
	return s;
}

}
}
