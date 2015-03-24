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
#include "fon.h"

namespace pagmo { namespace problem {

/**
 * Will construct Fonseca and Fleming's study problem.
 *
 * @see problem::base constructors.
 */
fon::fon():base(3,0,2)
{
	set_bounds(-4,4);
}

/// Clone method.
base_ptr fon::clone() const
{
	return base_ptr(new fon(*this));
}

/// Implementation of the objective function.
void fon::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 2 && x.size() == 3);
	f[0] = 1 - std::exp(-(x[0] - 1 / std::sqrt(3.0)) * (x[0] - 1 / std::sqrt(3.0)) - (x[1] - 1 / std::sqrt(3.0)) * (x[1] -
		1 / std::sqrt(3.0)) - (x[2] - 1 / std::sqrt(3.0)) * (x[2] - 1 / std::sqrt(3.0)));
	f[1] = 1 - std::exp(-(x[0] + 1 / std::sqrt(3.0)) * (x[0] + 1 / std::sqrt(3.0)) - (x[1] + 1 / std::sqrt(3.0)) * (x[1] +
		1 / std::sqrt(3.0)) - (x[2] + 1 / std::sqrt(3.0)) * (x[2] + 1 / std::sqrt(3.0)));
}

std::string fon::get_name() const
{
	return "Fonseca and Fleming's study";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::fon)
