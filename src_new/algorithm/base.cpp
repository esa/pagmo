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
#include <typeinfo>

#include "../rng.h"
#include "base.h"

namespace pagmo
{
namespace algorithm {

/// Default constructor.
/**
 * Will initialise the internal random number generators using seeds provided by the thread-safe pagmo::static_rng_uint32.
 */
base::base():m_drng(static_rng_uint32()()),m_urng(static_rng_uint32()()) {}

/// Trivial destructor.
base::~base() {}

/// Return human readable representation of the algorithm.
/**
 * Will return a formatted string containing the problem type (in mangled C++ form).
 * The output of human_readable_extra() will be appended at the end of the string.
 */
std::string base::human_readable() const
{
	std::ostringstream s;
	s << "Algorithm type: " << typeid(*this).name() << '\n';
	s << human_readable_extra();
	return s.str();
}

/// Extra information in human readable format.
/**
 * Default implementation returns an empty string.
 */
std::string base::human_readable_extra() const
{
	return std::string();
}

/// Overload stream operator for algorithm::base.
/**
 * Equivalent to printing base::human_readable() to stream.
 */
std::ostream &operator<<(std::ostream &s, const base &alg)
{
	s << alg.human_readable();
	return s;
}

}
}
