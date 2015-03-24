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

#include <cmath>
#include <boost/math/constants/constants.hpp>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "dtlz.h"

static int __check__(int N)
{
		if (N > 7 || N < 1) {
				pagmo_throw(value_error, "the problem id needs to be one of [1..7]");
		}
		return N;
}

namespace pagmo { namespace problem {

static const double PI_HALF = boost::math::constants::pi<double>() / 2.0;

/**
 * Will construct a problem of the DTLZ test-suite.
 * @param[in] id problem number: 1-6 possible.
 * @param[in] k paramter defining integer dimension of the problem: k + fdim - 1
 * @param[in] fdim number of objectives
 * @param[in] alpha controls density of solutions (used only by DTLZ4)
 *
 * @see problem::base constructors.
 */
dtlz::dtlz(size_type id, size_type k, size_type fdim, const size_t alpha)
	:base_dtlz(k + fdim - 1, fdim),
	m_problem_number(__check__(id)),
	m_alpha(alpha)
{
	// Set bounds.
	set_lb(0.0);
	set_ub(1.0);
}

/// Clone method.
base_ptr dtlz::clone() const
{
	return base_ptr(new dtlz(*this));
}

/// Convergence metric for a decision_vector (0 = converged to the optimal front)
double dtlz::g_func(const decision_vector &x) const
{
	switch(m_problem_number)
	{ // We start with the 6-7 cases as for absurd reasons behind my comprehension this is way more efficient
	case 6:
		return g6_func(x);
	case 7:
		return g7_func(x);
	case 1:
	case 3:
		return g13_func(x);
	case 2:
	case 4:
	case 5:
		return g245_func(x);
	default:
		pagmo_throw(value_error, "Error: There are only 7 test functions in this test suite!");
	}
	return -1.0;
}

/// Implementation of the objective function.
void dtlz::objfun_impl(fitness_vector &f, const decision_vector &x) const
{

	pagmo_assert(f.size() == get_f_dimension());
	pagmo_assert(x.size() == get_dimension());

	switch(m_problem_number)
	{
	case 1:
	f1_objfun_impl(f,x);
		break;
	case 2:
	case 3:
		f23_objfun_impl(f,x);
		break;
	case 4:
		f4_objfun_impl(f,x);
		break;
	case 5:
	case 6:
		f56_objfun_impl(f,x);
		break;
	case 7:
		f7_objfun_impl(f,x);
		break;
	default:
		pagmo_throw(value_error, "Error: There are only 7 test functions in this test suite!");
		break;
	}
}

/// Implementations of the different g-functions used
double dtlz::g13_func(const decision_vector &x) const
{
	double y = 0.0;
	for(decision_vector::size_type i = 0; i < x.size(); ++i) {
		y += pow(x[i] - 0.5, 2) - cos(20 * boost::math::constants::pi<double>() * (x[i] - 0.5));
	}
	return 100.0 * (y + x.size());
}

double dtlz::g245_func(const decision_vector &x) const
{
	double y = 0.0;
	for(decision_vector::size_type i = 0; i < x.size(); ++i) {
		y += pow(x[i] - 0.5, 2);
	}
	return y;
}

double dtlz::g6_func(const decision_vector &x) const
{
	double y = 0.0;
	for(decision_vector::size_type i = 0; i < x.size(); ++i) {
		y += pow(x[i], 0.1);
	}
	return y;
}

double dtlz::g7_func(const decision_vector &x) const
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
double dtlz::h7_func(const fitness_vector &f, const double g) const
{
	// NOTE: we intentionally ignore the last element of the vector to make things easier
	int fdim = get_f_dimension();
	double y = 0.0;

	for(decision_vector::size_type i = 0; i < f.size() - 1; ++i) {
		y += (f[i] / (1.0 + g)) * (1.0 + sin(3 * boost::math::constants::pi<double>() * f[i]) );
	}
	return fdim - y;
}

/// Implementation of the objective functions.
/* The chomosome: x_1, x_2, ........, x_M-1, x_M, .........., x_M+k
 *											 [------- Vector x_M -------]
 *               x[0], x[1], ... ,x[fdim-2], x[fdim-1], ... , x[fdim+k-1] */
void dtlz::f1_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
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

void dtlz::f23_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	// computing distance-function
	decision_vector x_M;
	double g;

	for(problem::base::size_type i = f.size() - 1; i < x.size(); ++i) {
		x_M.push_back(x[i]);
	}

	g = g_func(x_M);

	// computing shape-functions
	f[0] = (1.0 + g);
	for(problem::base::size_type i = 0; i < f.size() - 1; ++i) {
		f[0] *= cos(x[i] * PI_HALF);
	}

	for(problem::base::size_type i = 1; i < f.size() - 1; ++i) {
		f[i] = (1.0 + g);
		for(problem::base::size_type j = 0; j < f.size() - (i+1); ++j) {
			f[i] *= cos(x[j] * PI_HALF);
		}
		f[i] *= sin(x[f.size() - (i+1)] * PI_HALF);
	}

	f[f.size()-1] = (1.0 + g) * sin(x[0] * PI_HALF);

}

void dtlz::f4_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	// computing distance-function
	decision_vector x_M;
	double g;

	for(problem::base::size_type i = f.size() - 1; i < x.size(); ++i) {
		x_M.push_back(x[i]);
	}

	g = g_func(x_M);

	// computing shape-functions
	f[0] = (1.0 + g);
	for(problem::base::size_type i = 0; i < f.size() - 1; ++i) {
		f[0] *= cos(pow(x[i],m_alpha) * PI_HALF);
	}

	for(problem::base::size_type i = 1; i < f.size() - 1; ++i) {
		f[i] = (1.0 + g);
		for(problem::base::size_type j = 0; j < f.size() - (i+1); ++j) {
			f[i] *= cos(pow(x[j], m_alpha) * PI_HALF);
		}
		f[i] *= sin(pow(x[f.size() - (i+1)], m_alpha) * PI_HALF);
	}

	f[f.size()-1] = (1.0 + g) * sin(pow(x[0],m_alpha) * PI_HALF);

}

void dtlz::f56_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	// computing distance-function
	decision_vector x_M;
	double g;

	for(problem::base::size_type i = f.size() - 1; i < x.size(); ++i) {
		x_M.push_back(x[i]);
	}

	g = g_func(x_M);

	// computing meta-variables
	decision_vector theta(f.size(), 0.0);
	double t;

	theta[0] = x[0]; // * PI_HALF;
	t =  1.0 / (2.0 * (1.0 + g));

	for(problem::base::size_type i = 1; i < f.size(); ++i) {
		theta[i] = t + ((g * x[i]) / (1.0 + g));
	}

	// computing shape-functions
	f[0] = (1.0 + g);
	for(problem::base::size_type i = 0; i < f.size() - 1; ++i) {
		f[0] *= cos(theta[i] * PI_HALF);
	}

	for(problem::base::size_type i = 1; i < f.size() - 1; ++i) {
		f[i] = (1.0 + g);
		for(problem::base::size_type j = 0; j < f.size() - (i+1); ++j) {
			f[i] *= cos(theta[j] * PI_HALF);
		}
		f[i] *= sin(theta[f.size() - (i+1)] * PI_HALF);
	}

	f[f.size()-1] = (1.0 + g) * sin(theta[0] * PI_HALF);

}

void dtlz::f7_objfun_impl(fitness_vector &f, const decision_vector &x) const
{
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

	f[f.size()-1] = (1.0 + g) * h7_func(f, g);

}

std::string dtlz::get_name() const
{
	std::string retval("DTLZ");
	retval.append(boost::lexical_cast<std::string>(m_problem_number));

	return retval;
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::dtlz)
