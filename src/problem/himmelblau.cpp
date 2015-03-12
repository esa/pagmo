/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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
himmelblau::himmelblau():base(-6.,6.,2)
{
	// initialize best solution
	initialize_best();
}

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

/// Initialize the Himmelblau's function four local minima.
void himmelblau::initialize_best(void)
{
	const int min_count = 4;
	const int x_dimension = 2;

	const double x_vector[][x_dimension] = {
		{3.0, 2.0},
		{-2.805118, 3.131312},
		{-3.779310, -3.283186},
		{3.584428, -1.848126}
	};

	std::vector<decision_vector> best_x(min_count);

	for (int i = 0; i < min_count; ++i) {
		decision_vector x(x_dimension);
		std::copy(x_vector[i],x_vector[i] + x_dimension,x.begin());
		best_x[i] = x;
	}

	set_best_x(best_x);
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::himmelblau)
