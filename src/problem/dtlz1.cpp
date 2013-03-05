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
#include "../population.h"
#include "base.h"
#include "dtlz1.h"

namespace pagmo { namespace problem {

/**
 * Will construct DTLZ1.
 *
 * @param[in] k paramter defining integer dimension of the problem: k + fdim - 1
 * @param[in] fdim number of objectives
 *
 * @see problem::base constructors.
 */
dtlz1::dtlz1(int k, fitness_vector::size_type fdim):base(k + fdim - 1, 0, fdim)
{
	// Set bounds.
	set_lb(0.0);
	set_ub(1.0);
}

/// Clone method.
base_ptr dtlz1::clone() const
{
	return base_ptr(new dtlz1(*this));
}

/// Gives a convergence metric for the population (0 = converged to the optimal front)
double dtlz1::p_distance(const pagmo::population &pop) const
{
	double c = 0.0;
	f_size_type fdim = pop.problem().get_f_dimension();
	decision_vector x_M;
	decision_vector x;
	
    for (std::vector<double>::size_type i = 0; i < pop.size(); ++i) {
		x_M.clear();
		x = pop.get_individual(i).cur_x;
		for(problem::base::size_type j = fdim - 1; j < x.size(); ++j) {
			x_M.push_back(x[j]);
		}
		c += g_func(x_M);
    }

    return c / pop.size();
}

/// Implementation of the distance function g
double dtlz1::g_func(const decision_vector &x) const
{
	double y = 0.0;
	for(decision_vector::size_type i = 0; i < x.size(); ++i) {
		y += pow(x[i] - 0.5, 2) - cos(20 * boost::math::constants::pi<double>() * (x[i] - 0.5));
	}
	return 100.0 * (y + x.size());
}

/// Implementation of the objective function.
/* The chomosome: x_1, x_2, ........, x_M-1, x_M, .........., x_M+k
 *											 [------- Vector x_M -------]
 *               x[0], x[1], ... ,x[fdim-2], x[fdim-1], ... , x[fdim+k-1] */
void dtlz1::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == get_f_dimension());
    pagmo_assert(x.size() == get_dimension());

	// computing distance-function
	decision_vector x_M;
	double g;

	for(problem::base::size_type i = f.size() - 1; i < x.size(); ++i) {
		x_M.push_back(x[i]);
	}

	g = g_func(x_M);

	// computing shape-functions
	f[0] = 0.5 * (1.0 + g);
	
	for(problem::base::size_type i = 0; i < f.size() - 1; ++i) {
		f[0] *= x[i];
	}

	for(problem::base::size_type i = 1; i < f.size() - 1; ++i) {
		f[i] = 0.5 * (1.0 + g);
		for(problem::base::size_type j = 0; j < f.size() - (i+1); ++j) {
			f[i] *= x[j];
		}
		f[i] *= 1 - x[f.size() - (i+1)];
	}

	f[f.size()-1] = 0.5 *  (1 - x[0]) * (1.0 + g);

}

std::string dtlz1::get_name() const
{
	return "DTLZ1";
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::dtlz1);
