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
#include "rastrigin.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional Rastrigin problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
rastrigin::rastrigin(int n):base(n)
{
	// Set bounds.
	set_lb(-5.12);
	set_ub(5.12);
}

/// Clone method.
base_ptr rastrigin::clone() const
{
	return base_ptr(new rastrigin(*this));
}

/// Implementation of the objective function.
void rastrigin::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	const double omega = 2.0 * boost::math::constants::pi<double>();
	f[0] = 0;
	const decision_vector::size_type n = x.size();
	for (decision_vector::size_type i = 0; i < n; ++i) {
		f[0] += x[i] * x[i] - 10.0 * std::cos(omega * x[i]);
	}
	f[0] += 10.0 * n;
}

std::string rastrigin::get_name() const
{
	return "Rastrigin";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::rastrigin)
