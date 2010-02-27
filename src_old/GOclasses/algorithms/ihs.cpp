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

// 09/01/2009: Initial version by Francesco Biscani.

#include <cmath>
#include <iostream>
#include <vector>

#include "../../exceptions.h"
#include "../basic/population.h"
#include "../problems/base.h"
#include "base.h"
#include "ihs.h"

namespace pagmo
{
namespace algorithm {

/// Simple constructor.
/**
 * Initialises the parameters of the algorithms to "standard" values: phmcr = 0.85, ppar_min = 0.35, ppar_max = 0.99,
 * bw_min = 1E-5 and bw_max = 1.
 * @param[in] gen number of generations.
 * @throws value_error if gen is not positive.
 */
ihs::ihs(int gen):base(),m_gen(gen),m_phmcr(0.85),m_ppar_min(0.35),m_ppar_max(0.99),m_bw_min(1E-5),m_bw_max(1)
{
	if (gen <= 0) {
		pagmo_throw(value_error,"number of generations must be positive");
	}
}

/// Advanced constructor.
/**
 * Allows to specify in detail the parameters of the algorithm.
 * @param[in] gen number of generations.
 * @param[in] phmcr rate of choosing from memory.
 * @param[in] ppar_min minimum pitch adjustment rate.
 * @param[in] ppar_max maximum pitch adjustment rate.
 * @param[in] bw_min minimum distance bandwidth.
 * @param[in] bw_max maximum distance bandwidth.
 * @throws value_error if gen is not positive, phmcr is not in the ]0,1[ interval, ppar min/max are not in the ]0,1[ interval,
 * min/max quantities are less than/greater than max/min quantities.
 */
ihs::ihs(int gen, const double &phmcr, const double &ppar_min, const double &ppar_max,
         const double &bw_min, const double &bw_max):
		base(),m_gen(gen),m_phmcr(phmcr),m_ppar_min(ppar_min),m_ppar_max(ppar_max),m_bw_min(bw_min),m_bw_max(bw_max)
{
	if (gen <= 0) {
		pagmo_throw(value_error,"number of generations must be positive");
	}
	if (phmcr >=1 || phmcr <= 0 || ppar_min >=1 || ppar_min <= 0 || ppar_max >=1 || ppar_max <= 0) {
		pagmo_throw(value_error,"probabilities must be in the ]0,1[ range");
	}
	if (ppar_min >= ppar_max) {
		pagmo_throw(value_error,"minimum pitch adjustment rate must be smaller than maximum pitch adjustment rate");
	}
	if (bw_min <= 0 || bw_max <= bw_min) {
		pagmo_throw(value_error,"bandwidth values must be positive, and minimum bandwidth must be smaller than maximum bandwidth");
	}
}

/// Evolve population.
/**
 * @throws value_error if population is empty.
 */
population ihs::evolve(const population &pop) const
{
	// Let's store some useful variables.
	const problem::base &problem = pop.problem();
	const std::vector<double> &lb = problem.get_lb(), &ub = problem.get_ub();
	const size_t problem_size = lb.size(), pop_size = pop.size();
	if (pop_size == 0) {
		pagmo_throw(value_error,"cannot evolve an empty population");
	}
	// This is the return population.
	population retval = pop;
	// Temporary decision vector, and lower-upper bounds difference vector.
	std::vector<double> tmp_dv(problem_size), lu_diff(problem_size);
	for (size_t i = 0; i < problem_size; ++i) {
		lu_diff[i] = ub[i] - lb[i];
	}
	const double c = std::log(m_bw_min/m_bw_max) / m_gen;
	for (size_t g = 0; g < m_gen; ++g) {
		const double ppar_cur = m_ppar_min + ((m_ppar_max - m_ppar_min) * g) / m_gen, bw_cur = m_bw_max * std::exp(c * g);
		for (size_t i = 0; i < problem_size; ++i) {
			const double next_rn = drng();
			if (drng() <= m_phmcr) {
				// With random probability, tmp's i-th chromosome element is the one from a randomly chosen individual.
				tmp_dv[i] = retval[(size_t)(next_rn * pop_size)].get_decision_vector()[i];
				if (drng() <= ppar_cur) {
					// Randomly, add or subtract pitch from the current chromosome element.
					const double next_next_rn = drng();
					if (drng() > .5) {
						tmp_dv[i] += next_next_rn * bw_cur * lu_diff[i];
					} else {
						tmp_dv[i] -= next_next_rn * bw_cur * lu_diff[i];
					}
					// Handle the case in which we addded or subtracted too much and ended up out
					// of boundaries.
					if (tmp_dv[i] > ub[i]) {
						tmp_dv[i] = ub[i];
					} else if (tmp_dv[i] < lb[i]) {
						tmp_dv[i] = lb[i];
					}
				}
			} else {
				tmp_dv[i] = lb[i] + next_rn * lu_diff[i];
			}
		}
		const double tmp_fitness = problem.objfun(tmp_dv);
		const individual &worst = retval.extractWorstIndividual();
		if (tmp_fitness < worst.get_fitness()) {
			retval.replace_worst(individual(tmp_dv, worst.get_velocity(), tmp_fitness));
		}
	}
	return retval;
}

void ihs::log(std::ostream &s) const
{
	s << "IHS - generations:" << m_gen << " phmcr:" << m_phmcr << " ppar_min:" << m_ppar_min
	<< " ppar_max:" << m_ppar_max << " bw_min:" << m_bw_min << " bw_max:" << m_bw_max;
}

}
}
