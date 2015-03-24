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

#include "gtoc_1.h"
#include "../AstroToolbox/mga.h"

namespace pagmo { namespace problem {

/// Problem Constructor
/**
 * @see problem::base constructors.
 */
gtoc_1::gtoc_1():base(8),Delta_V(8),rp(6),t(8)
{
	// Set bounds.
	const double lb[8] = {3000,14,14,14,14,100,366,300};
	const double ub[8] = {10000,2000,2000,2000,2000,9000,9000,9000};
	set_bounds(lb,lb+8,ub,ub+8);

	//Defining the problem data up the problem parameters
	problem.type = asteroid_impact;
	problem.mass = 1500.0;				// Satellite initial mass [Kg]
	problem.Isp = 2500.0;               // Satellite specific impulse [s]
	problem.DVlaunch = 2.5;				// Launcher DV in km/s

	int sequence_[8] = {3,2,3,2,3,5,6,10}; // sequence of planets
	std::vector<int> sequence(8);
	problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_+8);

	const int rev_[8] = {0,0,0,0,0,0,1,0}; // sequence of clockwise legs
	std::vector<int> rev(8);
	problem.rev_flag.insert(problem.rev_flag.begin(), rev_, rev_+8);

	problem.asteroid.keplerian[0] = 2.5897261;    // Asteroid data
	problem.asteroid.keplerian[1] = 0.2734625;
	problem.asteroid.keplerian[2] = 6.40734;
	problem.asteroid.keplerian[3] = 128.34711;
	problem.asteroid.keplerian[4] = 264.78691;
	problem.asteroid.keplerian[5] = 320.479555;
	problem.asteroid.epoch = 53600;
}

/// Clone method.
base_ptr gtoc_1::clone() const
{
	return base_ptr(new gtoc_1(*this));
}

/// Implementation of the objective function.
void gtoc_1::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	MGA(x,problem,rp,Delta_V,f[0]);
}

std::string gtoc_1::get_name() const
{
	return "GTOC_1";
}



}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::gtoc_1)
