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

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <string>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "schwefel.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional Schwefel problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
schwefel::schwefel(int n):base(n)
{
	// Set bounds.
	set_lb(-500);
	set_ub(500);
}

/// Clone method.
base_ptr schwefel::clone() const
{
	return base_ptr(new schwefel(*this));
}

/// Implementation of the objective function.
void schwefel::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	std::vector<double>::size_type n = x.size();
	double value=0;

	for (std::vector<double>::size_type i=0; i<n; i++){
		value += x[i] * sin(sqrt(fabs(x[i])));
		}
		f[0] = 418.9828872724338 * n - value;
}

std::string schwefel::get_name() const
{
	return "Schwefel";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::schwefel)
