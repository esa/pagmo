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
#include "pressure_vessel.h"

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
 * Will construct the pressure vessel problem
 *
 */
pressure_vessel::pressure_vessel():base(4,0,1,4,4,__constraint_tolerances__(4,4))
{
	// initialize best solution
	initialize_best();

	// set the bounds for the current problem
	const double lb[] = {1.,1.,10.,10.};
	const double ub[] = {99.,99.,200.,200.};
	set_bounds(lb,ub);
}

/// Clone method.
base_ptr pressure_vessel::clone() const
{
	return base_ptr(new pressure_vessel(*this));
}

/// Implementation of the objective function.
void pressure_vessel::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	/* objective function */
	f[0] = 0.6224 * x[0] * x[2] * x[3] +
			1.7781 * x[1] * x[2]*x[2] +
			3.1661 * x[0]*x[0] * x[3] +
			19.84 * x[0]*x[0] * x[2];
}

/// Implementation of the constraint function.
void pressure_vessel::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	/* constraints g<=0 */
	c[0] = - x[0] + 0.0193 * x[2];
	c[1] = - x[1] + 0.00954 * x[2];
	c[2] = - boost::math::constants::pi<double>() * x[2]*x[2] * x[3] -
			(4./3.) * boost::math::constants::pi<double>() * x[2]*x[2]*x[2] + 1296000.;
	c[3] = x[3] - 240.;
}

std::string pressure_vessel::get_name() const
{
	std::string retval("Pressure vessel");
	return retval;
}

void pressure_vessel::initialize_best(void)
{
	std::vector<decision_vector> best_x;

	int x_dimension = 4;
	// Coello and Montes
	const double x_vector[] = {0.812500, 0.437500, 42.097398, 176.654050};

	decision_vector x(x_dimension);
	std::copy(x_vector,x_vector + x_dimension,x.begin());
	best_x.push_back(x);

	set_best_x(best_x);
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::pressure_vessel)
