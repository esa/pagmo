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

#include <vector>

#include "../../exceptions.h"
#include "../basic/population.h"
#include "../problems/GOproblem.h"
#include "hs_algorithm.h"

hs_algorithm::hs_algorithm(int gen, const double &phmcr, const double &ppar, const double &bw):
	m_gen(gen),m_phmcr(phmcr),m_ppar(ppar),m_bw(bw)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (phmcr >=1 || phmcr <= 0 || ppar >= 1 || ppar <= 0 || bw >= 1 || bw <= 0) {
		pagmo_throw(value_error,"probabilities and bandwidth must be in the ]0,1[ range");
	}
}

Population hs_algorithm::evolve(const Population &pop) const
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
	// Temporary decision vector.
	std::vector<double> tmp_dv(problem_size);
	// Vector of picks.
	std::vector<size_t> picks;
	picks.reserve(problem_size);
	for (size_t g = 0; g < m_gen; ++g) {
		for (size_t i = 0; i < problem_size; ++i) {
			const double next_rn = drng();
			if (drng() <= m_phmcr) {
				const size_t pick = (size_t)(next_rn * pop_size);
				picks.push_back(pick);
				tmp_dv[i] = pop[pick].getDecisionVector()[i];
			} else {
				tmp_dv[i] = lb[i] + next_rn * (ub[i] - lb[i]);
			}
		}
		const size_t picks_size = picks.size();
		for (size_t i = 0; i < picks_size; ++i) {
			if (drng() <= m_ppar) {
				const double next_rn = drng();
				if (drng() > .5) {
					tmp_dv[picks[i]] += next_rn * m_bw * (ub[picks[i]] - lb[picks[i]]);
				} else {
					tmp_dv[picks[i]] -= next_rn * m_bw * (ub[picks[i]] - lb[picks[i]]);
				}
				if (tmp_dv[picks[i]] > ub[picks[i]]) {
					tmp_dv[picks[i]] = ub[picks[i]];
				} else if (tmp_dv[picks[i]] < lb[picks[i]]) {
					tmp_dv[picks[i]] = lb[picks[i]];
				}
			}
		}
		const double tmp_fitness = problem.objfun(tmp_dv);
		Individual &worst = retval.worst();
		if (tmp_fitness < worst.getFitness()) {
			worst = Individual(tmp_dv,tmp_fitness);
		}
	}
	return retval;
}
