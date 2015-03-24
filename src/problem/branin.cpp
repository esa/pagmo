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

#include "../types.h"
#include "base.h"
#include "branin.h"

namespace pagmo { namespace problem {

/// Constructor.
branin::branin():base(2,0,1,0,0)
{
	const double lb[] = {-5,0};
	const double ub[] = {10,15};
	set_bounds(lb,ub);
}

/// Clone method.
base_ptr branin::clone() const
{
	return base_ptr(new branin(*this));
}

/// Implementation of the objective function.
void branin::objfun_impl(fitness_vector &fv, const decision_vector &x) const
{
	const double x1 = x[0], x2 = x[1];
	const double a = 1, b = 5.1 / (4 * boost::math::constants::pi<double>() * boost::math::constants::pi<double>()),
		c = 5 / boost::math::constants::pi<double>(), d = 6, e = 10, f = 1 / (8 * boost::math::constants::pi<double>());
	fv[0] = a * (x2 - b * x1 * x1 + c * x1 - d) * (x2 - b * x1 * x1 + c * x1 - d) + e * (1 - f) * std::cos(x1) + e;
}

std::string branin::get_name() const
{
	return "Branin";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::branin)
