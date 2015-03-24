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

#include <algorithm>
#include <vector>
#include <set>

#include "../population.h"
#include "base.h"
#include "base_s_policy.h"
#include "hv_greedy_s_policy.h"
#include "best_s_policy.h"
#include "../exceptions.h"
#include "../util/hypervolume.h"

using namespace pagmo::util;

namespace pagmo { namespace migration {

/// Constructor from migration rate and type.
/**
 * @param[in] rate migration rate.
 * @param[in] type migration rate type.
 * @param[in] nadir_eps epsilon value for the nadir point computation
 *
 * @see base_s_policy::base_s_policy.
 */
hv_greedy_s_policy::hv_greedy_s_policy(const double &rate, rate_type type, const double nadir_eps):base_s_policy(rate,type), m_nadir_eps(nadir_eps) { }

base_s_policy_ptr hv_greedy_s_policy::clone() const
{
	return base_s_policy_ptr(new hv_greedy_s_policy(*this));
}

std::vector<population::individual_type> hv_greedy_s_policy::select(population &pop) const
{
	// Fall back to best_s_policy when facing a single-objective problem.
	if (pop.problem().get_f_dimension() == 1) {
		return best_s_policy(m_rate, m_type).select(pop);
	}

	pagmo_assert(get_n_individuals(pop) <= pop.size());
	// Gets the number of individuals to select
	const population::size_type migration_rate = get_n_individuals(pop);
	// Create a temporary array of individuals.
	std::vector<population::individual_type> result;

	// Indices of fronts.
	std::vector< std::vector< population::size_type> > fronts_i = pop.compute_pareto_fronts();

	// Fitness vectors of individuals according to the indices above.
	std::vector< std::vector< fitness_vector> > fronts_f (fronts_i.size());

	// Nadir point is established manually later, first point is as a first "safe" candidate.
	fitness_vector refpoint(pop.get_individual(0).cur_f);

	for (unsigned int f_idx = 0 ; f_idx < fronts_i.size() ; ++f_idx) {
		fronts_f[f_idx].resize(fronts_i[f_idx].size());
		for (unsigned int p_idx = 0 ; p_idx < fronts_i[f_idx].size() ; ++p_idx) {
			fronts_f[f_idx][p_idx] = fitness_vector(pop.get_individual(fronts_i[f_idx][p_idx]).cur_f);

			// Update the nadir point manually for efficiency.
			for (unsigned int d_idx = 0 ; d_idx < fronts_f[f_idx][p_idx].size() ; ++d_idx) {
				refpoint[d_idx] = std::max(refpoint[d_idx], fronts_f[f_idx][p_idx][d_idx]);
			}
		}
	}

	// Epsilon is added to nadir point
	for (unsigned int d_idx = 0 ; d_idx < refpoint.size() ; ++d_idx) {
		refpoint[d_idx] += m_nadir_eps;
	}

	// Store which front we process (start with front 0) and the number of processed individuals.
	unsigned int front_idx = 0;
	unsigned int processed_individuals = 0;

	// Vector for maintaining the original indices of points
	std::vector<unsigned int> orig_indices;

	while (processed_individuals < migration_rate) {
		// If we need to pull every point from given front anyway, just push back the individuals right away
		if (fronts_f[front_idx].size() <= (migration_rate - processed_individuals)) {
			for(unsigned int i = 0 ; i < fronts_i[front_idx].size() ; ++i) {
				result.push_back(pop.get_individual(fronts_i[front_idx][i]));
			}

			processed_individuals += fronts_f[front_idx].size();
			++front_idx;
		} else {
			// Prepare the vector for the original indices
			if (orig_indices.size() == 0) {
				orig_indices.resize(fronts_i[front_idx].size());
				iota(orig_indices.begin(), orig_indices.end(), 0);
			}

			// Compute the greatest contributor
			hypervolume hv(fronts_f[front_idx], false);
			hv.set_copy_points(false);
			unsigned int gc_idx = hv.greatest_contributor(refpoint);
			result.push_back(pop.get_individual(fronts_i[front_idx][orig_indices[gc_idx]]));
			
			// Remove it from the front along with its index
			orig_indices.erase(orig_indices.begin() + gc_idx);
			fronts_f[front_idx].erase(fronts_f[front_idx].begin() + gc_idx);
			++processed_individuals;
		}
	}

	return result;
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::migration::hv_greedy_s_policy)
