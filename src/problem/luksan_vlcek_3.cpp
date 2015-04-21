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
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "luksan_vlcek_3.h"

static int __check__(int N){
	if (N < 8 || (N) % 4)
	{
		pagmo_throw(value_error,"problem dimension N needs to be at least 8 and a multiple of 4");
	}
	return N;
}

namespace pagmo { namespace problem {

/// Constructor.
/**
 * Construct the problem from its dimension.
 * Setting cub=clb=0 creates an instance of the original Luksan Vlcek equality constrained problem, example  5.3
 * Using clb < cub allows the obtain a problem formulation with inequality constraints
 *
 * @param[in] N Problem dimension
 * @param[in] clb lower bounds for the constraints.
 * @param[in] cub upper bounds for the constraints.
 * @throws value_error if N is smaller than 6 and N+2 (the resulting problem dimension)
 * is not multiple of 4, cub < clb
 *
 * @see L.Luksan and J.Vlcek, "Sparse and Parially Separable Test Problems for Unconstrained and Equality Constrained Optimization"
 */
luksan_vlcek_3::luksan_vlcek_3(int N, const double &clb, const double &cub):base(__check__(N), 0, 1, 2 * 2, 2 * 2, 1e-14)
{
	if (clb > cub)
	{
		pagmo_throw(value_error,"constraints lower bound is higher than the upper bound");
	}
	set_lb(-5);
	set_ub(5);
	m_clb = std::vector<double>(2,clb);
	m_cub = std::vector<double>(2,cub);
}

/// Clone method.
base_ptr luksan_vlcek_3::clone() const
{
	return base_ptr(new luksan_vlcek_3(*this));
}

/// Implementation of the objective function.
void luksan_vlcek_3::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = 0.;
	for (decision_vector::size_type i=0; i<(x.size()-2)/2; i++)
	{
		double a1 = x[2*i]+10.*x[2*i+1];
		double a2 = x[2*i+2] - x[2*i+3];
		double a3 = x[2*i+1] - 2.*x[2*i+2];
		double a4 = x[2*i] - x[2*i+3];
		f[0] += a1*a1 + 5.*a2*a2 + std::pow(a3,4)+ 10.*std::pow(a4,4);
	}

}
/// Implementation of the constraint function.
void luksan_vlcek_3::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	int n = x.size();
	c[0] = 3.*std::pow(x[0],3) + 2.*x[1] - 5. + std::sin(x[0]-x[1])*std::sin(x[0]+x[1]) - m_cub[0];
	c[1] = m_clb[0] - ( 3.*std::pow(x[0],3) + 2.*x[1] - 5. + std::sin(x[0]-x[1])*std::sin(x[0]+x[1]) );
	c[2] = 4.*x[n-3] - x[n-4]*std::exp(x[n-4]-x[n-3]) - 3 - m_cub[1];
	c[3] = m_clb[1] - ( 4.*x[n-3] - x[n-4]*std::exp(x[n-4]-x[n-3]) - 3 );
}

/// Implementation of the sparsity structure
void luksan_vlcek_3::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	//Initial point
	decision_vector x0(get_dimension(),1);
	//Numerical procedure
	estimate_sparsity(x0, lenG, iGfun, jGvar);
}

std::string luksan_vlcek_3::get_name() const
{
	return "Luksan-Vlcek 3";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::luksan_vlcek_3)
