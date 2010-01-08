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

// 09/09/09 Created by Francesco Biscani.

#include <boost/lexical_cast.hpp>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "../../exceptions.h"
#include "../../Functions/rng/rng.h"
#include "../basic/population.h"
#include "base.h"
#include "sp_test.h"

namespace pagmo
{
namespace problem {

inventory::inventory(int sample_size):base(1),d(sample_size),m_sample_size(sample_size)
{
	set_lb(0,0);
	set_ub(0,100);
}


double inventory::objfun_(const std::vector<double> &x) const
{
	const double c=1.0,b=1.5,h=0.1;
	double retval=0;
	for (size_t i = 0; i<m_sample_size; ++i) {
		retval += c * x[0] + b * std::max(d[i]-x[0],0.0) + h * std::max(x[0] - d[i], 0.0);
	}
	return retval / m_sample_size;
}

void inventory::pre_evolution(population & pop) const
{
	for (size_t i = 0; i<m_sample_size; ++i) {
		d[i] = static_rng_double()() * 100;
	}
	//Re-evaluate the population with respect to the new seed (Internal Sampling Method)
	for (size_t i=0; i<pop.size(); ++i) {
		pop[i] = individual(*this, pop[i].get_decision_vector(), pop[i].get_velocity());
	}
}

}
}
