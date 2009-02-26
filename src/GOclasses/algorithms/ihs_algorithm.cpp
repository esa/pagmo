/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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
#include <vector>

#include "../../exceptions.h"
#include "../basic/population.h"
#include "../problems/GOproblem.h"
#include "ihs_algorithm.h"

ihs_algorithm::ihs_algorithm(int gen, const double &phmcr, const double &ppar_min, const double &ppar_max,
	const double &bw_min, const double &bw_max):
	m_gen(gen),m_phmcr(phmcr),m_ppar_min(ppar_min),m_ppar_max(ppar_max),m_bw_min(bw_min),m_bw_max(bw_max)
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

Population ihs_algorithm::evolve(const Population &pop) const
{
	// Let's store some useful variables.
	const GOProblem &problem = pop.problem();
	const std::vector<double> &lb = problem.getLB(), &ub = problem.getUB();
	const size_t problem_size = lb.size(), pop_size = pop.size();
	if (pop_size == 0) {
		pagmo_throw(index_error,"cannot evolve an empty population");
	}
	// This is the return population.
	Population retval = pop;
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
				tmp_dv[i] = retval[(size_t)(next_rn * pop_size)].getDecisionVector()[i];
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
		const Individual &worst = retval.extractWorstIndividual();
		if (tmp_fitness < worst.getFitness()) {
			retval.replace_worst(Individual(tmp_dv,worst.getVelocity(),tmp_fitness));
		}
	}
	return retval;
}

void ihs_algorithm::log(std::ostream &s) const
{
	s << "IHS - generations:" << m_gen << " phmcr:" << m_phmcr << " ppar_min:" << m_ppar_min
		<< " ppar_max:" << m_ppar_max << " bw_min:" << m_bw_min << " bw_max:" << m_bw_max;
}
