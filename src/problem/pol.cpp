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

#include <boost/math/constants/constants.hpp>
#include <cmath>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "pol.h"

namespace pagmo { namespace problem {

/**
 * Will construct Poloni's study problem.
 *
 * @see problem::base constructors.
 */
pol::pol():base(2,0,2)
{
	// Set bounds.
	set_lb(-std::atan(1.0f) * 4.0f); //-PI
	set_ub(std::atan(1.0f) * 4.0f); //PI
}

/// Clone method.
base_ptr pol::clone() const
{
	return base_ptr(new pol(*this));
}

/// Implementation of the objective function.
void pol::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 2);
	pagmo_assert(x.size() == 2);

	f[0] = 0;
	f[1] = 0;

	double a1 = 0.5 * sin(1.0) - 2 * cos(1.0) + sin(2.0) - 1.5 * cos(2.0);
 
	double a2 = 1.5 * sin(1.0) - cos(1.0) + 2 * sin(2.0) - 0.5 * cos(2.0);
 
        double b1 = 0.5 * sin(x[0]) - 2 * cos(x[0]) + sin(x[1]) - 1.5 * cos(x[1]);
 
  	double b2 = 1.5 * sin(x[0]) - cos(x[0]) + 2 * sin(x[1]) - 0.5 * cos(x[1]);

	f[0] = 1 + (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2);
	f[1] = (x[0]+3)*(x[0]+3) + (x[1]+1)*(x[1]+1);
}

std::string pol::get_name() const
{
	return "Poloni's study";
}


}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::pol)
