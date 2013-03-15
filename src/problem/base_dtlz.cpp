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

#include "base_dtlz.h"
#include "../types.h"
#include "../population.h"

namespace pagmo { namespace problem {

/// Constructor
base_dtlz::base_dtlz(int n, int nf):base(n, 0, nf, 0, 0, 0.0) {
	};

/// Gives a convergence metric for the population (0 = converged to the optimal front)
double base_dtlz::p_distance(const pagmo::population &pop) const
{
	double c = 0.0;
	f_size_type fdim = pop.problem().get_f_dimension();
	decision_vector x_M;
	decision_vector x;
	
    for (std::vector<double>::size_type i = 0; i < pop.size(); ++i) {
		x_M.clear();
		x = pop.get_individual(i).cur_x;
		for(problem::base::size_type j = fdim - 1; j < x.size(); ++j) {
			x_M.push_back(x[j]);
		}
		c += g_func(x_M);
    }

    return c / pop.size();
}

/// Implementation of the distance function g
double base_dtlz::g_func(const decision_vector &x) const
{
	(void) x;	// to avoid warnings during compilation, parameter x is just used by the base classes
	std::cout << "No g-function implemented, return 0.0 as function value" << std::endl;
	return 0.0;
}

}} //namespaces
