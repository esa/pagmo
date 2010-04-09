/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

#include "messenger_full.h"
#include "../AstroToolbox/mga_dsm.h"

namespace pagmo { namespace problem {

/// Problem Constructor
/**
 * @see problem::base constructors.
 */
const int sequence[7] = {3, 2, 2, 1, 1, 1, 1};
messenger_full::messenger_full():base(26),problem(orbit_insertion,sequence,7,0,0,0,0.704,2440 + 200)
{
	// Set bounds.
	const double lb[26] = {1900, 3,    0, 0, 100, 100, 100, 100, 100, 100, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI};
	const double ub[26] = {2200, 4.05, 1, 1, 500, 500, 500, 500, 500, 550, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,   6,   6,    6,    6,    6,  M_PI,  M_PI,  M_PI,  M_PI,  M_PI};
	set_bounds(lb,lb+26,ub,ub+26);

	// Set sequence

}

/// Clone method.
base_ptr messenger_full::clone() const
{
	return base_ptr(new messenger_full(*this));
}

/// Implementation of the objective function.
void messenger_full::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	MGA_DSM(x, problem,f[0]);
}

}}
