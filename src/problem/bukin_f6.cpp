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
#include "bukin_f6.h"

namespace pagmo { namespace problem {

/// Default constructor.
bukin::bukin():base(2)
{
	// Set bounds.
	set_lb(0,-15.0);
	set_ub(0,-5.0);
	set_lb(1,-3.0);
	set_ub(1,3.0);
}

/// Clone method.
base_ptr bukin::clone() const
{
	return base_ptr(new bukin(*this));
}

/// Implementation of the objective function.
void bukin::objfun_impl(fitness_vector &f, const decision_vector &vx) const
{
	pagmo_assert(f.size() == 1 && vx.size() == 2);
	const double x = vx[0], y = vx[1];
	f[0] = 100 * sqrt(fabs(y - 0.01 * x * x)) + 0.01*fabs(x+10);
}

std::string bukin::get_name() const
{
	return "Bukin f6";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::bukin)
