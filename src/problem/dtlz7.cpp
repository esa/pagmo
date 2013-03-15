/*****************************************************************************
 *   Copyright (C) 2004-2012 The PaGMO development team,                     *
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

#include <cmath>
#include <boost/math/constants/constants.hpp>

#include "../exceptions.h"
#include "../types.h"
#include "base_dtlz.h"
#include "dtlz7.h"

namespace pagmo { namespace problem {

/**
 * Will construct DTLZ7.
 *
 * @param[in] k paramter defining integer dimension of the problem: k + fdim - 1
 * @param[in] fdim number of objectives
 *
 * @see problem::base constructors.
 */
dtlz7::dtlz7(int k, fitness_vector::size_type fdim):base_dtlz(k + fdim - 1, fdim)
{
	// Set bounds.
	set_lb(0.0);
	set_ub(1.0);
}

/// Clone method.
base_ptr dtlz7::clone() const
{
	return base_ptr(new dtlz7(*this));
}

/// Implementation of the distance function g
double dtlz7::g_func(const decision_vector &x) const
{
	// NOTE: the original g-function should return 1 + (9.0 / x.size()) * y but we drop the 1
	// to have the minimum at 0.0 so we can use the p_distance implementation in base_dtlz
	// to have the p_distance converging towards 0.0 rather then towards 1.0
	double y = 0.0;
	for(decision_vector::size_type i = 0; i < x.size(); ++i) {
		y += x[i];
	}
	return (9.0 / x.size()) * y;
}

/// Implementation of the distribution function h
double dtlz7::h_func(const fitness_vector &f, const double g) const
{
	// NOTE: we intentionally ignore the last element of the vector to make things easier
	int fdim = get_f_dimension();
	double y = 0.0;

	for(decision_vector::size_type i = 0; i < f.size() - 1; ++i) {
		y += (f[i] / (1.0 + g)) * (1.0 + sin(3 * boost::math::constants::pi<double>() * f[i]) );
	}
	return fdim - y;
}


/// Implementation of the objective function.
/* The chomosome: x_1, x_2, ........, x_M-1, x_M, .........., x_M+k
 *											 [------- Vector x_M -------]
 *               x[0], x[1], ... ,x[fdim-2], x[fdim-1], ... , x[fdim+k-1] */
void dtlz7::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == get_f_dimension());
    pagmo_assert(x.size() == get_dimension());

	// computing distance-function
	decision_vector x_M;
	double g;

	for(problem::base::size_type i = f.size() - 1; i < x.size(); ++i) {
		x_M.push_back(x[i]);
	}

	g = 1.0 + g_func(x_M);		// +1.0 according to the original definition of the g-function for DTLZ7

	// computing shape-functions
	for(problem::base::size_type i = 0; i < f.size() - 1; ++i) {
		f[i] = x[i];
	}

	f[f.size()-1] = (1.0 + g) * h_func(f, g);

}

std::string dtlz7::get_name() const
{
	return "DTLZ7";
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::dtlz7);
