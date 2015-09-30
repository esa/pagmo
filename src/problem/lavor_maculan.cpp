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

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include <cmath>
#include "lavor_maculan.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an lavor_maculan problem with hydrocarbon chain of N atoms (N-3 parameters).
 *
 * @param[in] atoms number of atoms
 *
 * @see problem::base constructors.
 */
lavor_maculan::lavor_maculan(int atoms): base(atoms - 3)
{
	if (atoms <= 0 || atoms < 4) {
		pagmo_throw(value_error, "number of atoms for lavor-maculan problem must be positive and greater than 3");
	}
	// Set bounds.
	set_lb(0.0);
	set_ub(5.0);

	// Initialise global minima vector
	std::vector<decision_vector> best_x(1);
	// Single global minimum of Lavor Maculan is a repeating sequence of
	// (1.039195303, 3.141592654, 1.039195303, 3.141592654, 1.039195303, ...)
	double rep[] = {1.039195303, 3.141592654};
	for(int i = 0; i < atoms - 3; ++i) {
		best_x[0].push_back(rep[i % 2]);
	}

	set_best_x(best_x);
}

/// Clone method.
base_ptr lavor_maculan::clone() const
{
	return base_ptr(new lavor_maculan(*this));
}

/// Implementation of the objective function.
void lavor_maculan::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	const decision_vector::size_type n = x.size();
	f[0] = 0.0;
	for (decision_vector::size_type i = 0 ; i < n ; ++i) {
		f[0] += 1 + cos(3 * x[i]) + (((i % 2 == 1) ? 1 : -1) / (sqrt(10.60099896 - 4.141720682 * cos(x[i]))));
	}
}

std::string lavor_maculan::get_name() const
{
	return "Lavor Maculan";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::lavor_maculan)
