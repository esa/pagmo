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
#include "../exceptions.h"
#include "base_unc_mo.h"

namespace pagmo { namespace problem {

/// Constructor from dimension and fitness dimension
/**
 * Will construct an n dimensional unconstrained multi-objective problem
 * with nf objectives
 *
 * @param[in] n dimension of the problem.
 * @param[in] ni integer dimension of the problem.
 * @param[in] nf number of objectives
 *
 * @see problem::base constructors.
 */
base_unc_mo::base_unc_mo(size_type n, size_type ni, f_size_type nf):base(n, ni, nf, 0, 0, 0.0) {}

/// Distance from the Pareto front (of a population)
/**
 * Will return the average across the entire population of the convergence metric
 *
 * @param[in] pop population to be assigned a pareto distance
 *
 * @see problem::base_unc_mo::p_distance virtual method.
 */
double base_unc_mo::p_distance(const pagmo::population &pop) const
{
	double c = 0.0;
	for (population::size_type i = 0; i < pop.size(); ++i) {
		c += convergence_metric(pop.get_individual(i).cur_x);
	}

	return c / pop.size();
}

/// Distance from the Pareto front (of a decision_vector)
/**
 * Will return the convergence metric of the decision_vector
 *
 * @param[in] x decision_vector 
 *
 */
double base_unc_mo::p_distance(const decision_vector &x) const
{
	return convergence_metric(x);
}

/// Default implementation for a convergence metric
/**
 *
 * @param[in] x decision_vector 
 *
 * @throws not_implemented_error always
 *
 */
double base_unc_mo::convergence_metric(const decision_vector &x) const
{
	(void) x;	// avoids warnings during compilation
	pagmo_throw(not_implemented_error, "Error: a convergence metric is not implemented for this problem.");
}

}} // namespaces
