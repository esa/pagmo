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
#include "himmelblau.h"

namespace pagmo { namespace problem {

/// Default constructor.
himmelblau::himmelblau():base(-6.,6.,2) {}

/// Clone method.
base_ptr himmelblau::clone() const
{
	return base_ptr(new himmelblau(*this));
}

/// Implementation of the objective function.
void himmelblau::objfun_impl(fitness_vector &f, const decision_vector &xv) const
{
	pagmo_assert(f.size() == 1 && xv.size() == get_dimension());
	const double x = xv[0], y = xv[1];
	f[0] = ((x * x + y - 11) * (x * x + y - 11) + (x + y * y - 7) * (x + y * y - 7));
}

std::string himmelblau::get_name() const
{
	return "Himmelblau";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::himmelblau)
