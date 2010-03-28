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
#include <cmath>
#include <string>

#include "../types.h"
#include "base.h"
#include "string_match_mo.h"

namespace pagmo { namespace problem {

/// Constructor from C++ string.
string_match_mo::string_match_mo(const std::string &str):base(boost::numeric_cast<int>(str.size()),boost::numeric_cast<int>(str.size()),
	boost::numeric_cast<problem::base::f_size_type>(str.size()),0,0),m_str(str)
{
	set_bounds(0,255);
}

/// Constructor from C string.
string_match_mo::string_match_mo(const char *str):base(boost::numeric_cast<int>(std::string(str).size()),boost::numeric_cast<int>(std::string(str).size()),
	boost::numeric_cast<problem::base::f_size_type>(std::string(str).size()),0,0),m_str(str)
{
	set_bounds(0,255);
}

base_ptr string_match_mo::clone() const
{
	return base_ptr(new string_match_mo(*this));
}

void string_match_mo::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	for (std::string::size_type i = 0; i < m_str.size(); ++i) {
		f[boost::numeric_cast<problem::base::f_size_type>(i)] =
			std::abs(static_cast<char>(x[boost::numeric_cast<size_type>(i)]) - m_str[i]);
	}
}

std::string string_match_mo::human_readable_extra() const
{
	return std::string("\tString: \"") + m_str + "\"";
}

} }
