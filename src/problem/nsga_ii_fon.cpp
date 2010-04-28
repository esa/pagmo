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
#include <string>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "nsga_ii_fon.h"

namespace pagmo { namespace problem {

nsga_ii_fon::nsga_ii_fon():base(3,0,2,0,0)
{
	set_bounds(-4,4);
}

base_ptr nsga_ii_fon::clone() const
{
	return base_ptr(new nsga_ii_fon(*this));
}

void nsga_ii_fon::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 2 && x.size() == 3);
	f[0] = 1 - std::exp(-(x[0] - 1 / std::sqrt(3)) * (x[0] - 1 / std::sqrt(3)) - (x[1] - 1 / std::sqrt(3)) * (x[1] -
		1 / std::sqrt(3)) - (x[2] - 1 / std::sqrt(3)) * (x[2] - 1 / std::sqrt(3)));
	f[1] = 1 - std::exp(-(x[0] + 1 / std::sqrt(3)) * (x[0] + 1 / std::sqrt(3)) - (x[1] + 1 / std::sqrt(3)) * (x[1] +
		1 / std::sqrt(3)) - (x[2] + 1 / std::sqrt(3)) * (x[2] + 1 / std::sqrt(3)));
}

std::string nsga_ii_fon::get_name() const
{
	return "NSGA-II FON";
}

}}
