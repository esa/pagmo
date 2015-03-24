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
base::base():m_screen_output(false), m_drng(rng_generator::get<rng_double>()), m_urng(rng_generator::get<rng_uint32>()), m_fevals(0) {}

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

/// Sets screen output
/**
 * Sets screen output boolean variable. When True the algorithm may print stuff on screen (careful, in multithreading this can mess up things)
 *
 * @param[in] p true or false
 */
void base::set_screen_output(const bool p) {m_screen_output = p;}

/// Gets screen output
/**
 * Gets screen output boolean variable. When True the algorithm may print stuff on screen 
 * (careful, in multithreading this can mess up things)
 *
 * return boolean value associated to the screen output option
 */
bool base::get_screen_output() const {return m_screen_output;}

/// Resets the random number generators
/**
 * Resets the seed of the internal random number generators to p
 *
 * @param[in] p the new seed. It will be used to instantiate both drng and urng 
 */
void base::reset_rngs(const unsigned int p) const {
	m_urng = rng_uint32(p);
	m_drng = rng_double(p);
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
	s << "Algorithm name: " << get_name();
	const std::string tmp(human_readable_extra());
	if (tmp.size()) {
		s << " - " << tmp;
	}
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
