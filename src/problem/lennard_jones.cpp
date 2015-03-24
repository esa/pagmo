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
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "lennard_jones.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct a Lennard-Jones problem
 *
 * @param[in] atoms number of atoms
 *
 * @see problem::base constructors.
 */
lennard_jones::lennard_jones(int atoms):base(3*atoms-6)
{
	if (atoms <= 0 || atoms < 3) {
		pagmo_throw(value_error,"number of atoms for lennard-jones problem must be positive and greater than 2");
	}
	for (int i = 0; i < 3*atoms-6; i++) {
		if ( (i != 0) && (i % 3) == 0 ) {
			set_lb(i,0.0);
			set_ub(i,6.0);
		} else {
			set_lb(i,-3.0);
			set_ub(i,3.0);
		}
	}
}

/// Clone method.
base_ptr lennard_jones::clone() const
{
	return base_ptr(new lennard_jones(*this));
}

/// Helper function that transforms the decision vector x in atoms positions r
double lennard_jones::r(const int& atom, const int& coord, const std::vector<double>& x) {
	if(atom == 0) { //x1,y1,z1 fixed
		return 0.0;
	} else if(atom == 1) {
		if(coord < 2) { //x2,y2    fixed
			return 0.0;
		} else { //z2 is a variable
			return x[0];
		}
	} else if(atom == 2) {
		if(coord == 0) { //x3	   fixed
			return 0.0;
		} else { //y3 and x3 are variables
			return x[coord];
		}
	} else {
		return x[3 * (atom - 2) + coord];
	}
}

/// Implementation of the objective function.
void lennard_jones::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	std::vector<double>::size_type n = x.size();
	int atoms = (n + 6) / 3;
	double sixth, dist;

	f[0] = 0;
	//We evaluate the potential
	for ( int i=0; i<(atoms-1); i++ ) {
		for ( int j=(i+1); j<atoms; j++ ) {
			dist = pow(r(i, 0, x) - r(j, 0, x), 2) + pow(r(i, 1, x) - r(j, 1, x), 2) + pow(r(i, 2, x) - r(j, 2, x), 2);  //rij^2
			if ( dist == 0.0 ) {
				f[0] = 1e+20;	//penalty
			}
			else {
				sixth = pow(dist, -3);	//rij^-6
				f[0] += (pow(sixth, 2) - sixth);
			}
		}
	}
	f[0] = 4 * f[0];
}

std::string lennard_jones::get_name() const
{
	return "Lennard-Jones";
}

}}//namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::lennard_jones)
