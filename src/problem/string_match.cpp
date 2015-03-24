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

#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <string>

#include "../types.h"
#include "base.h"
#include "string_match.h"

namespace pagmo { namespace problem {

/// Constructor from C++ string.
string_match::string_match(const std::string &str):base(boost::numeric_cast<int>(str.size()),boost::numeric_cast<int>(str.size()),1,0,0),m_str(str)
{
	set_bounds(0,255);
}

/// Constructor from C string.
string_match::string_match(const char *str):base(boost::numeric_cast<int>(std::string(str).size()),boost::numeric_cast<int>(std::string(str).size()),1,0,0),m_str(str)
{
	set_bounds(0,255);
}

base_ptr string_match::clone() const
{
	return base_ptr(new string_match(*this));
}

void string_match::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	double retval = 0;
	for (std::string::size_type i = 0; i < m_str.size(); ++i) {
		retval += std::abs(static_cast<char>(x[boost::numeric_cast<size_type>(i)]) - m_str[i]);
	}
	f[0] = retval;
}


/// Returns the encoded string
/**
 * @return the encoded string (std::string)
 */
std::string string_match::pretty(const decision_vector &x) const {
	std::string retval;
	for (decision_vector::size_type i = 0; i < x.size(); ++i) {
		retval += char(x[i]);
	}
	return retval;
}

std::string string_match::human_readable_extra() const
{
	return std::string("\tString: \"") + m_str + "\"\n";
}

std::string string_match::get_name() const
{
	return "String match";
}

} }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::string_match)
