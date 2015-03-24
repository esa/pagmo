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

#include "sagas.h"
#include "../AstroToolbox/mga_dsm.h"

namespace pagmo { namespace problem {

const int sagas::sequence[3] = {3, 3, 5};

/// Constructor
/**
* Instantiates the rosetta problem*/
sagas::sagas():base(12), problem(time2AUs,sequence,3,0,0,0,0,0)
{
	//Setting the objective parameters
	problem.AUdist = 50.0;
	problem.DVtotal = 6.782;
	problem.DVonboard = 1.782;

	const double lb[12] = {7000, 0, 0, 0, 50, 300, 0.01, 0.01, 1.05, 8, -M_PI, -M_PI};
	const double ub[12] = {9100, 7, 1, 1, 2000, 2000, 0.9, 0.9 ,7 ,500 ,M_PI,M_PI};
	set_bounds(lb,lb+12,ub,ub+12);
}

/// Clone method.
base_ptr sagas::clone() const
{
	return base_ptr(new sagas(*this));
}

/// Implementation of the objective function.
void sagas::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	MGA_DSM(x, problem,f[0]);
}

/// Implementation of the sparsity structure.
/**
 * No sparsity present (box-constrained problem).
 */

void sagas::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	lenG=12;
	iGfun.resize(lenG);
	jGvar.resize(lenG);
	for (int i = 0; i<lenG; ++i)
	{
		iGfun[i] = 0;
		jGvar[i] = i;
	}
}

std::string sagas::get_name() const
{
	return "Sagas";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::sagas)
