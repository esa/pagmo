/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

#include "base.h"
#include "luksan_vlcek_1.h"

namespace pagmo { namespace problem {

luksan_vlcek_1::luksan_vlcek_1(size_t N, const double clb, const double cub):base(N,0,1,2*(N-2),2*(N-2))
{
	if (N <=2)
	{
		pagmo_throw(value_error,"Problem dimension needs to be at least 3");
	}
	if (clb >cub)
	{
		pagmo_throw(value_error,"constraints lower bound is higher than the upper bound");
	}
	set_lb(-5);
	set_ub(5);
	m_clb = std::vector<double>(N-2,clb);
	m_cub = std::vector<double>(N-2,cub);
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
		c[2 * i] =  (3.*pow(x[i+1],3.) + 2.*x[i+2] - 5.
		+ sin(x[i+1]-x[i+2])*sin(x[i+1]+x[i+2]) + 4.*x[i+1]
		- x[i]*exp(x[i]-x[i+1]) - 3.) - m_cub[i];
		c[2 * i + 1] = - (3.*pow(x[i+1],3.) + 2.*x[i+2] - 5.
		+ sin(x[i+1]-x[i+2])*sin(x[i+1]+x[i+2]) + 4.*x[i+1]
		- x[i]*exp(x[i]-x[i+1]) - 3.) + m_clb[i];
	}
}

/// Implementation of the sparsity structure: automated detection
void luksan_vlcek_1::set_sparsity(int& lenG, std::vector<int>& iGfun, std::vector<int>& jGvar) const
{
	//Initial point
	decision_vector x0(get_dimension(),1);
	//Numerical procedure
	estimate_sparsity(x0, lenG, iGfun, jGvar);
}

} }
