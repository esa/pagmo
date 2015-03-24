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

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "tens_comp_string.h"

static const std::vector<double> __constraint_tolerances__(int c_dimension, int ic_dimension)
{
	std::vector<double> constraint_tolerances(c_dimension);
	// equality constraints
	for(int i=0; i<c_dimension-ic_dimension; i++) {
		constraint_tolerances[i] = 0.0001;
	}
	// inequality constraints
	for(int i=c_dimension-ic_dimension; i<c_dimension; i++) {
		constraint_tolerances[i] = 0.;
	}
	return constraint_tolerances;
}

namespace pagmo { namespace problem {

/// Constructor
/**
 * Will construct the tension compression string problem
 *
 */
tens_comp_string::tens_comp_string():base(3,0,1,4,4,__constraint_tolerances__(4,4))
{
	// initialize best solution
	initialize_best();

	// set the bounds for the current problem
	const double lb[] = {0.05,0.25,2.};
	const double ub[] = {2.,1.3,15.};
	set_bounds(lb,ub);
}

/// Clone method.
base_ptr tens_comp_string::clone() const
{
	return base_ptr(new tens_comp_string(*this));
}

/// Implementation of the objective function.
void tens_comp_string::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	/* objective function */
	f[0] = x[0]*x[0] * x[1] * (x[2] + 2.);
}

/// Implementation of the constraint function.
void tens_comp_string::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	c[0] = 1. - (x[1]*x[1]*x[1] * x[2]) / (71785. * x[0]*x[0]*x[0]*x[0]);
	c[1] = (4. * x[1]*x[1] - x[0] * x[1]) /
			( 12566. * x[0]*x[0]*x[0] * (x[1] - x[0]) ) +
			1./(5108. * x[0]*x[0]) - 1.;
	c[2] = 1. - 140.45 * x[0] / (x[2] * x[1]*x[1]);
	c[3] = (x[0] + x[1]) / 1.5 - 1.;
}

std::string tens_comp_string::get_name() const
{
	std::string retval("Tension compression string");
	return retval;
}

void tens_comp_string::initialize_best(void)
{
	std::vector<decision_vector> best_x;

	int x_dimension = 3;
	// He and Wang
	const double x_vector[] = {0.051728, 0.357644, 11.244543};

	decision_vector x(x_dimension);
	std::copy(x_vector,x_vector + x_dimension,x.begin());
	best_x.push_back(x);

	set_best_x(best_x);
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::tens_comp_string)
