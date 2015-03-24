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

#include "messenger.h"
#include "../AstroToolbox/mga_dsm.h"
#include "../AstroToolbox/misc4Tandem.h"

namespace pagmo { namespace problem {

const int messenger::sequence[5] = {3,3,2,2,1};

/// Constructor
/**
* Instantiates the rosetta problem*/
messenger::messenger():base(18), problem(total_DV_rndv,sequence,5,0,0,0,0,0)
{

	const double lb[18] = {1000, 1, 0, 0, 200, 30,  30,  30,  0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.1, -M_PI, -M_PI, -M_PI};
	const double ub[18] = {4000, 5, 1, 1, 400, 400, 400, 400, 0.99, 0.99, 0.99, 0.99, 6,   6,   6,    M_PI,  M_PI,  M_PI};
	set_bounds(lb,lb+18,ub,ub+18);
}

/// Clone method.
base_ptr messenger::clone() const
{
	return base_ptr(new messenger(*this));
}

/// Implementation of the objective function.
void messenger::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	MGA_DSM(x, problem,f[0]);
}

/// Implementation of the sparsity structure.
/**
 * No sparsity present (box-constrained problem).
 */

void messenger::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	lenG=18;
	iGfun.resize(lenG);
	jGvar.resize(lenG);
	for (int i = 0; i<lenG; ++i)
	{
		iGfun[i] = 0;
		jGvar[i] = i;
	}
}

std::string messenger::get_name() const
{
	return "Messenger";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::messenger)
