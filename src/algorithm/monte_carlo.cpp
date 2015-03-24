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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <cstddef>

#include "../population.h"
#include "../types.h"
#include "base.h"
#include "monte_carlo.h"

namespace pagmo { namespace algorithm {

/// Constructor from number of iterations.
monte_carlo::monte_carlo(int n):base(),m_max_eval(boost::numeric_cast<std::size_t>(n)) {}

/// Clone method.
base_ptr monte_carlo::clone() const
{
	return base_ptr(new monte_carlo(*this));
}

/// Evolve method.
void monte_carlo::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type prob_dimension = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type pop_size = pop.size();
	// Get out if there is nothing to do.
	if (pop_size == 0 || m_max_eval == 0) {
		return;
	}
	// Initialise temporary decision vector, fitness vector and decision vector.
	decision_vector tmp_x(prob_dimension);
	fitness_vector tmp_f(prob.get_f_dimension());
	constraint_vector tmp_c(prob.get_c_dimension());
	// Main loop.
	for (std::size_t i = 0; i < m_max_eval; ++i) {
		// Generate a random decision vector.
		for (problem::base::size_type j = 0; j < prob_dimension - prob_i_dimension; ++j) {
			tmp_x[j] = boost::uniform_real<double>(lb[j],ub[j])(m_drng);
		}
		for (problem::base::size_type j = prob_dimension - prob_i_dimension; j < prob_dimension; ++j) {
			tmp_x[j] = boost::uniform_int<int>(lb[j],ub[j])(m_urng);
		}
		// Compute fitness and constraints.
		prob.objfun(tmp_f,tmp_x);
		prob.compute_constraints(tmp_c,tmp_x);
		// Locate the worst individual.
		const population::size_type worst_idx = pop.get_worst_idx();
		if (prob.compare_fc(tmp_f,tmp_c,pop.get_individual(worst_idx).cur_f,pop.get_individual(worst_idx).cur_c)) {
			pop.set_x(worst_idx,tmp_x);
		}
	}
}

/// Algorithm name
std::string monte_carlo::get_name() const
{
	return "Monte Carlo";
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string monte_carlo::human_readable_extra() const
{
	std::ostringstream s;
	s << "max_eval:" << m_max_eval;
	return s.str();
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::monte_carlo)
