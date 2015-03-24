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

#include <string>

#include "base.h"
#include "snopt_toyprob.h"

namespace pagmo { namespace problem {

/// Default constructor.
/**
 * Search bounds are set to \f$ x \in \left[ 0,10 \right] \f$ and \f$ y \in \left[ -10,10 \right]\f$.
 */
snopt_toyprob::snopt_toyprob():base(2,0,1,2,2)
{
	const double lb[] = {0,-10};
	const double ub[] = {10,10};
	set_lb(lb);
	set_ub(ub);
}

/// Clone method.
base_ptr snopt_toyprob::clone() const
{
	return base_ptr(new snopt_toyprob(*this));
}

/// Implementation of the objective function.
void snopt_toyprob::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = x[1];
}

/// Implementation of the constraint function.
void snopt_toyprob::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	c[0] = x[0] * x[0] + 4 * x[1] * x[1] - 4;
	c[1] = (x[0] - 2) * (x[0] - 2) + x[1] * x[1] - 5;
}

/// Implementation of the sparsity structure
void snopt_toyprob::set_sparsity(int& lenG, std::vector<int>& iGfun, std::vector<int>& jGvar) const
{
	//Initial point
	decision_vector x0(2);
	x0[0] = 1; x0[1] = 1;
	//Numerical procedure
	this->estimate_sparsity(x0, lenG, iGfun, jGvar);
}

std::string snopt_toyprob::get_name() const
{
	return "SNOPT toy problem";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::snopt_toyprob)
