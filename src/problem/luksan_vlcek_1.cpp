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

#include <boost/integer_traits.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "luksan_vlcek_1.h"

static int __check__(int N){
	if (N - 2 >= boost::integer_traits<int>::const_max / 2) {
		pagmo_throw(std::overflow_error,"overflow error");
	}
	if (N <= 2)
	{
		pagmo_throw(value_error,"problem dimension needs to be at least 3");
	}
	return N;
}

namespace pagmo { namespace problem {
/// Constructor.
/**
 * Construct the problem from its dimension.
 * Setting cub=clb=0 creates an instance of the original Luksan Vlcek equality constrained problem, example  5.1.
 * Using clb < cub allows the obtain a problem formulation with inequality constraints.
 *
 * @param[in] N problem dimension
 * @param[in] clb lower bounds for the constraints.
 * @param[in] cub upper bounds for the constraints.
 * @throws value_error if N is smaller than 2, cub < clb
 *
 * @see L. Luksan and J. Vlcek, "Sparse and Parially Separable Test Problems for Unconstrained and Equality Constrained Optimization"
 */
luksan_vlcek_1::luksan_vlcek_1(int N, const double &clb, const double &cub):base(__check__(N), 0, 1, 2 * (__check__(N) - 2), 2 * (__check__(N) - 2), 1e-14)
{

	if (clb > cub)
	{
		pagmo_throw(value_error,"constraints lower bound is higher than the upper bound");
	}
	set_lb(-5);
	set_ub(5);
	m_clb = decision_vector(boost::numeric_cast<decision_vector::size_type>(N-2),clb);
	m_cub = decision_vector(boost::numeric_cast<decision_vector::size_type>(N-2),cub);
}

/// Clone method.
base_ptr luksan_vlcek_1::clone() const
{
	return base_ptr(new luksan_vlcek_1(*this));
}

/// Implementation of the objective function.
void luksan_vlcek_1::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = 0.;
	for (pagmo::decision_vector::size_type i=0; i<x.size()-1; i++)
	{
		double a1 = x[i]*x[i]-x[i+1];
		double a2 = x[i] - 1.;
		f[0] += 100.*a1*a1 + a2*a2;
	}
}

/// Implementation of the constraint function.
void luksan_vlcek_1::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	for (pagmo::decision_vector::size_type i=0; i<x.size()-2; i++)
	{
		c[2 * i] =  (3.*std::pow(x[i+1],3.) + 2.*x[i+2] - 5.
		+ std::sin(x[i+1]-x[i+2])*std::sin(x[i+1]+x[i+2]) + 4.*x[i+1]
		- x[i]*std::exp(x[i]-x[i+1]) - 3.) - m_cub[i];
		c[2 * i + 1] = - (3.*std::pow(x[i+1],3.) + 2.*x[i+2] - 5.
		+ std::sin(x[i+1]-x[i+2])*std::sin(x[i+1]+x[i+2]) + 4.*x[i+1]
		- x[i]*std::exp(x[i]-x[i+1]) - 3.) + m_clb[i];
	}
}

/// Implementation of the sparsity structure: automated detection
void luksan_vlcek_1::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	//Initial point
	decision_vector x0(get_dimension(),1);
	//Numerical procedure
	estimate_sparsity(x0, lenG, iGfun, jGvar);
}

std::string luksan_vlcek_1::get_name() const
{
	return "Luksan-Vlcek 1";
}

} }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::luksan_vlcek_1)
