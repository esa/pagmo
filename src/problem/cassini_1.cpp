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

#include "cassini_1.h"
#include "../exceptions.h"
#include "../AstroToolbox/mga.h"

namespace pagmo { namespace problem {

/// Problem Constructor
/**
 * 
 * @param[in] objectives when equal to 1, the problem will be a single objective problem (DV)
 * When equal to 2, the problem is instantiated as a multiple-objectiove problem (DV,DT)
 * 
 * @see problem::base constructors.
 */
cassini_1::cassini_1(unsigned int objectives):base(6,0,objectives),Delta_V(6),rp(4),t(6)
{
	if (objectives != 1 && objectives !=2) {
		pagmo_throw(value_error,"Cassini_1 problem has either one or two objectives");
	}
	// Set bounds.
	const double lb[6] = {-1000,  30, 100, 30 , 400 , 1000};
	const double ub[6] = {0    , 400, 470, 400, 2000, 6000};
	set_bounds(lb,lb+6,ub,ub+6);

	// Set the mgaproblem data
	problem.type = total_DV_orbit_insertion; //Optimization type

	int sequence_[6] = {3,2,2,3,5,6}; //Sequence of planets
	problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_+6);

	const int rev_[6] = {0,0,0,0,0,0}; //Sequence of clockwise legs
	problem.rev_flag.insert(problem.rev_flag.begin(), rev_, rev_+6);

	problem.e  = 0.98;      // Insertion orbit eccentricity
	problem.rp = 108950;    // Insertion orbit pericenter
	problem.DVlaunch = 0;   // Launcher DV
}

/// Clone method.
base_ptr cassini_1::clone() const
{
	return base_ptr(new cassini_1(*this));
}

/// Implementation of the objective function.
void cassini_1::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	MGA(x,problem,rp,Delta_V,f[0]);
	if (get_f_dimension() == 2) {
		f[1] = (x[2]+x[3]+x[4]+x[5]); // + std::max(0.0,f[0] - 20) * 365.25;
	}	
	}

std::string cassini_1::get_name() const
{
	return "Cassini 1";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cassini_1)
