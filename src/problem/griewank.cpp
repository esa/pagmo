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

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "griewank.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional Griewank problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
griewank::griewank(int n):base(n)
{
	// Set bounds.
	set_lb(-600);
	set_ub(600);
}

/// Clone method.
base_ptr griewank::clone() const
{
	return base_ptr(new griewank(*this));
}

/// Implementation of the objective function.
void griewank::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	decision_vector::size_type n = x.size();
	double fr=4000.0;
	double retval = 0.0;
	double p = 1.0;

	for (decision_vector::size_type i=0; i<n; i++){ retval += x[i]*x[i];}
	for (decision_vector::size_type i=0; i<n; i++){ p *= cos(x[i]/sqrt(i+1.0));}
	f[0] = (retval/fr - p + 1);
}

std::string griewank::get_name() const
{
	return "Griewank";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::griewank)
