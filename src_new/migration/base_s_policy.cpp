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

#include <boost/numeric/conversion/cast.hpp>
#include <climits>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>

#include "../exceptions.h"
#include "../population.h"
#include "../types.h"
#include "base_s_policy.h"

namespace pagmo { namespace migration {

/// Constructor from migration rate and type.
/**
 * If migration type is absolute, then the input rate will be converted to the nearest integer. Will fail if
 * rate is outside the [0,INT_MAX] range.
 *
 * If migration type is fractional, the constructor will fail if rate is outside the [0,1] range.
 *
 * @param[in] rate migration rate.
 * @param[in] rate_type migration rate type.
 */
base_s_policy::base_s_policy(const double &rate,  migration_rate_type rate_type):m_rate(rate),m_rate_type(rate_type)
{
	switch (m_rate_type) {
		case absolute:
			if (m_rate < 0 || m_rate > INT_MAX) {
				pagmo_throw(value_error,"invalid absolute migration rate");
			}
			// Convert to nearest integer.
			m_rate = double_to_int::convert(m_rate);
			break;
		case fractional:
			if (m_rate < 0 || m_rate > 1) {
				pagmo_throw(value_error,"invalid fractional migration rate");
			}
			break;
		default:
			pagmo_throw(value_error,"invalid type for migration rate");
	}
}

/// Trivial destructor.
/**
 * No side effects.
 */
base_s_policy::~base_s_policy() {}

/// Get number of individuals to migrate from input population.
/**
 * @param[in] pop input population.
 *
 * @return the number of individuals to be migrated from the population according to the current migration rate and type.
 */
population::size_type base_s_policy::get_n_individuals(const population &pop) const
{
	population::size_type retval;
	switch (m_rate_type) {
		case absolute:
			retval = boost::numeric_cast<population::size_type>(m_rate);
			if (retval > pop.size()) {
				pagmo_throw(value_error,"absolute migration rate is higher than population size");
			}
			break;
		case fractional:
			retval = boost::numeric_cast<population::size_type>(double_to_int::convert(static_cast<double>(m_rate) * pop.size()));
			pagmo_assert(retval <= pop.size());
	}
	return retval;
}

/// Return human readable representation.
/**
 * Return a formatted string containing:
 * - the name of the policy (in C++ mangled form),
 * - migration type,
 * - migration rate,
 * - the output of human_readable_extra().
 *
 * @return human readable representation of the policy.
 */
std::string base_s_policy::human_readable() const
{
	std::ostringstream oss;
	oss << "Selection policy: " << typeid(*this).name() << '\n';
	oss << "\tMigration type: " << ((m_rate_type) ? "fractional" : "absolute") << '\n';
	oss << "\tMigration rate: " << ((m_rate_type) ? (m_rate * 100) : m_rate) << ((m_rate_type) ? "%" : "") << '\n';
	oss << human_readable_extra();
	return oss.str();
}

/// Return extra information in human readable representation.
/**
 * Return policy-specific information in human readable format. Default implementation returns an empty string.
 *
 * @return empty string.
 */
std::string base_s_policy::human_readable_extra() const
{
	return std::string();
}

/// Overload stream operator for migration::base_s_policy.
/**
 * Equivalent to printing base_s_policy::human_readable() to stream.
 *
 * @param[out] s std::ostream to which the policy will be streamed.
 * @param[in] sp migration::base_s_policy to be inserted into the stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const base_s_policy &sp)
{
	s << sp.human_readable();
	return s;
}

}}
