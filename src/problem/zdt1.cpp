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

#include <cmath>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "zdt1.h"

namespace pagmo { namespace problem {

/**
 * Will construct ZDT1.
 *
 * @see problem::base constructors.
 */
zdt1::zdt1():base(30,0,2)
{
	// Set bounds.
	set_lb(0.0);
	set_ub(1.0);
}

/// Clone method.
base_ptr zdt1::clone() const
{
	return base_ptr(new zdt1(*this));
}

/// Implementation of the objective function.
void zdt1::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 2);
	pagmo_assert(x.size() == 30);

	double g = 0;

	f[0] = x[0];

	for(problem::base::size_type i = 2; i < 30; ++i) {
		g += x[i];
	}
	g = 1 + (9 * g) / 29;
	
	f[1] = g * ( 1 - sqrt(x[0]/g));
	
}

std::string zdt1::get_name() const
{
	return "ZDT1";
}
}}
