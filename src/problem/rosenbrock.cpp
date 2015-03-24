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

#include <string>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "rosenbrock.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional Rosenbrock problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
rosenbrock::rosenbrock(int n):base(n)
{
	// Set bounds.
	set_lb(-5.0);
	set_ub(10);
}

/// Clone method.
base_ptr rosenbrock::clone() const
{
	return base_ptr(new rosenbrock(*this));
}

/// Implementation of the objective function.
void rosenbrock::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	const decision_vector::size_type n = x.size();
	f[0]=0;
	for (decision_vector::size_type i=0; i<n-1; ++i){
		f[0] += 100 * (x[i]*x[i] -x[i+1])*(x[i]*x[i] -x[i+1]) + (x[i]-1)*(x[i]-1);
	}
}

std::string rosenbrock::get_name() const
{
	return "Rosenbrock";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::rosenbrock)
