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

#include "../types.h"
#include "../population.h"
#include "base_dtlz.h"

namespace pagmo { namespace problem {

/// Constructor
base_dtlz::base_dtlz(int n, int nf):base_unc_mo(n, 0, nf) {}

/// Gives a convergence metric for the population (0 = converged to the optimal front)
double base_dtlz::convergence_metric(const decision_vector &x) const
{
	double c = 0.0;
	decision_vector x_M;
	for(problem::base::size_type j = get_f_dimension() - 1; j < x.size(); ++j) {
		x_M.push_back(x[j]);
	}
	c += g_func(x_M);
	return c;
}
}} //namespaces
