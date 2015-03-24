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
#include "welded_beam.h"

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
 * Will construct the welded beam problem
 *
 */
welded_beam::welded_beam():base(4,0,1,7,7,__constraint_tolerances__(7,7))
{
	// initialize best solution
	initialize_best();

	// set the bounds for the current problem
	const double lb[] = {0.1,0.1,0.1,0.1};
	const double ub[] = {2.,10.,10.,2.};
	set_bounds(lb,ub);
}

/// Clone method.
base_ptr welded_beam::clone() const
{
	return base_ptr(new welded_beam(*this));
}

/// Implementation of the objective function.
void welded_beam::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	/* objective function */
	f[0] = 1.10471 * x[0]*x[0] * x[1] + 0.04811 * x[2] * x[3] * (14. + x[1]);
}

/// Implementation of the constraint function.
void welded_beam::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	double P = 6000.;
	double L = 14.;
	double E = 30e+6;
	double G = 12e+6;
	double t_max = 13600.;
	double s_max = 30000.;
	double d_max = 0.25;

	double M = P*(L + x[1] / 2.);
	double R = std::sqrt(0.25 * (x[1]*x[1] + (x[0] + x[2]) * (x[0] + x[2])));
	double J = 2. / std::sqrt(2.) * x[0] * x[1] * (x[1]*x[1] / 12. + 0.25 *
			(x[0] + x[2])*(x[0] + x[2]));
	double P_c = (4.013 * E / (6. * L*L)) * x[2] * x[3]*x[3]*x[3] *
			(1 - 0.25 * x[2] * std::sqrt(E / G) / L);
	double t1 = P / (std::sqrt(2.) * x[0] * x[1]);
	double t2 = M * R / J;
	double t = sqrt(t1*t1 + t1 * t2 * x[1] / R + t2*t2);
	double s = 6. * P * L / (x[3] * x[2]*x[2]);
	double d = 4. * P * L*L*L / (E*x[3] * x[2]*x[2]*x[2]);

	c[0] = t - t_max;
	c[1] = s - s_max;
	c[2] = x[0] - x[3];
	c[3] = 0.10471 * x[0]*x[0] + 0.04811 * x[2] * x[3] * (14.0 + x[1]) - 5.0;
	c[4] = 0.125 - x[0];
	c[5] = d - d_max;
	c[6] = P - P_c;
}

std::string welded_beam::get_name() const
{
	std::string retval("Welded beam");
	return retval;
}

void welded_beam::initialize_best(void)
{
	std::vector<decision_vector> best_x;

	int x_dimension = 4;
	// Coello and Montes
	const double x_vector[] = {0.202369, 3.544214, 9.048210, 0.205723};

	decision_vector x(x_dimension);
	std::copy(x_vector,x_vector + x_dimension,x.begin());
	best_x.push_back(x);

	set_best_x(best_x);
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::welded_beam)
