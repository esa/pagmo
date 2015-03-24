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
#include "ackley.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional Schwefel problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
ackley::ackley(int n):base(n)
{
	// Set bounds.
	set_lb(-15.0);
	set_ub(30);

	// Initialise global minima vector vector vector vector
	std::vector<decision_vector> best_x(1);
	// The Ackley problem has only one global minimum: at (0, ..., 0)
	// (the value for all dimensions is 0).
	for(int i = 0; i < n; ++i) {
		best_x[0].push_back(0);
	}

	set_best_x(best_x);
}

/// Clone method.
base_ptr ackley::clone() const
{
	return base_ptr(new ackley(*this));
}

/// Implementation of the objective function.
void ackley::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	std::vector<double>::size_type n = x.size();
	f[0]=0;

	double omega = 2.0 * M_PI;
	double s1=0.0, s2=0.0;
	double nepero=exp(1.0);

	for (std::vector<double>::size_type i=0; i<n; i++){
		s1 += x[i]*x[i];
		s2 += cos(omega*x[i]);
	}
	f[0] = -20*exp(-0.2 * sqrt(1.0/n * s1))-exp(1.0/n*s2)+ 20 + nepero;
}

std::string ackley::get_name() const
{
	return "Ackley";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::ackley)
