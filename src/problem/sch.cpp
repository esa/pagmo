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
#include "sch.h"

namespace pagmo { namespace problem {

/**
 * Will construct Shaffer's study problem.
 *
 * @see problem::base constructors.
 */
sch::sch():base(1,0,2)
{
	set_bounds(-1000,1000);
}
/// Clone method.
base_ptr sch::clone() const
{
	return base_ptr(new sch(*this));
}

/// Implementation of the objective function.
void sch::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 2 && x.size() == 1);
	f[0] = x[0]*x[0];
	f[1] = (x[0]-2) * (x[0]-2);
}

std::string sch::get_name() const
{
	return "Shaffer's Study";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::sch)
