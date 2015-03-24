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

#include "rosetta.h"
#include "../AstroToolbox/mga_dsm.h"
#include "../AstroToolbox/misc4Tandem.h"

namespace pagmo { namespace problem {

const int rosetta::sequence[6] = {3, 3, 4, 3, 3, 10};

/// Constructor
/**
* Instantiates the rosetta problem*/
rosetta::rosetta():base(22), problem(rndv,sequence,6,0,0,0,0,0)
{
	//Setting the final comet elements
	problem.asteroid.keplerian[0] = 3.50294972836275;
	problem.asteroid.keplerian[1] = 0.6319356;
	problem.asteroid.keplerian[2] = 7.12723;
	problem.asteroid.keplerian[3] = 50.92302;
	problem.asteroid.keplerian[4] = 11.36788;
	problem.asteroid.keplerian[5] = 0.0;
	problem.asteroid.epoch = 52504.23754000012;
	problem.asteroid.mu = 0.0;

	const double lb[22] = {1460, 3, 0, 0, 300, 150, 150, 300, 700 , 0.01, 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI, -M_PI};
	const double ub[22] = {1825, 5, 1, 1, 500, 800, 800, 800, 1850, 0.9, 0.9 , 0.9 , 0.9 , 0.9 , 9   , 9   , 9   , 9   , M_PI , M_PI , M_PI , M_PI};
	set_bounds(lb,lb+22,ub,ub+22);
}

/// Clone method.
base_ptr rosetta::clone() const
{
	return base_ptr(new rosetta(*this));
}

/// Implementation of the objective function.
void rosetta::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	MGA_DSM(x, problem,f[0]);
}

/// Implementation of the sparsity structure.
/**
 * No sparsity present (box-constrained problem).
 */

void rosetta::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	lenG=22;
	iGfun.resize(22);
	jGvar.resize(22);
	for (int i = 0; i<lenG; ++i)
	{
		iGfun[i] = 0;
		jGvar[i] = i;
	}
}

std::string rosetta::get_name() const
{
	return "Rosetta";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::rosetta)
