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
#include <boost/numeric/conversion/cast.hpp>
#include <utility>
#include <vector>
#include <algorithm>

#include "../population.h"
#include "base.h"
#include "base_r_policy.h"
#include "hv_fair_r_policy.h"
#include "fair_r_policy.h"

#include "../util/hypervolume.h"


using namespace pagmo::util;

namespace pagmo { namespace migration {

/// Constructor from rate and rate type.
/**
 * @param[in] rate migration rate.
 * @param[in] type migration rate type.
 * @param[in] nadir_eps epsilon value for the nadir point computation
 *
 * @see base_r_policy::base_r_policy.
 */
hv_fair_r_policy::hv_fair_r_policy(const double &rate, rate_type type, const double nadir_eps):base_r_policy(rate,type), m_nadir_eps(nadir_eps) {}

base_r_policy_ptr hv_fair_r_policy::clone() const
{
	return base_r_policy_ptr(new hv_fair_r_policy(*this));
}

// Selection implementation.
std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> >
	hv_fair_r_policy::select(const std::vector<population::individual_type> &immigrants, const population &dest) const
{
	// Fall back to fair_r_policy when facing a single-objective problem.
	if (dest.problem().get_f_dimension() == 1) {
		return fair_r_policy(m_rate, m_type).select(immigrants, dest);
	}

	std::vector<population::individual_type> filtered_immigrants;
	filtered_immigrants.reserve(immigrants.size());

	// Keeps information on the original indexing of immigrants after we filter out the duplicates
	std::vector<unsigned int> original_immigrant_indices;
	original_immigrant_indices.reserve(immigrants.size());

	// Remove the duplicates from the set of immigrants
	std::vector<population::individual_type>::iterator im_it = (const_cast<std::vector<population::individual_type> &>(immigrants)).begin();
	unsigned int im_idx = 0;
	for( ; im_it != immigrants.end() ; ++im_it) {
		decision_vector im_x((*im_it).cur_x);

		bool equal = true;
		for ( unsigned int idx = 0 ; idx < dest.size() ; ++idx ) {
			decision_vector isl_x(dest.get_individual(idx).cur_x);
			equal = true;
			for (unsigned int d_idx = 0 ; d_idx < im_x.size() ; ++d_idx) {
				if (im_x[d_idx] != isl_x[d_idx]) {
					equal = false;
					break;
				}
			}
			if (equal) {
				break;
			}
		}
		if (!equal) {
			filtered_immigrants.push_back(*im_it);
			original_immigrant_indices.push_back(im_idx);
		}
		++im_idx;
	}

	// Computes the number of immigrants to be selected (accounting for the destination pop size)
	const population::size_type rate_limit = std::min<population::size_type>(get_n_individuals(dest), boost::numeric_cast<population::size_type>(filtered_immigrants.size()));

	// Defines the retvalue
	std::vector<std::pair<population::size_type, std::vector<population::individual_type>::size_type> > result;

	// Skip the remaining computation if there's nothing to do
	if (rate_limit == 0) {
		return result;
	}

	// Makes a copy of the destination population
	population pop_copy(dest);

	// Merge the immigrants to the copy of the destination population
	for (population::size_type i  = 0; i < rate_limit; ++i) {
		pop_copy.push_back(filtered_immigrants[i].cur_x);
	}

	// Population fronts stored as indices of individuals.
	std::vector< std::vector<population::size_type> > fronts_i = pop_copy.compute_pareto_fronts();

	// Population fronts stored as fitness vectors of individuals.
	std::vector< std::vector<fitness_vector> > fronts_f (fronts_i.size());

	// Nadir point is established manually later, first point is a first "safe" candidate.
	fitness_vector refpoint(pop_copy.get_individual(0).cur_f);

	// Fill fronts_f with fitness vectors and establish the nadir point
	for (unsigned int f_idx = 0 ; f_idx < fronts_i.size() ; ++f_idx) {
		fronts_f[f_idx].resize(fronts_i[f_idx].size());
		for (unsigned int p_idx = 0 ; p_idx < fronts_i[f_idx].size() ; ++p_idx) {
			fronts_f[f_idx][p_idx] = fitness_vector(pop_copy.get_individual(fronts_i[f_idx][p_idx]).cur_f);

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

	// Vector for maintaining the original indices of points for augmented population as 0 and 1
	std::vector<unsigned int> g_orig_indices(pop_copy.size(), 1);

	unsigned int no_discarded_immigrants = 0;

	// Store which front we process (start with the last front) and the number of processed individuals.
	unsigned int front_idx = fronts_i.size(); // front_idx is equal to the size, since it's decremented right in the main loop
	unsigned int processed_individuals = 0;

	// Pairs of (islander index, islander exclusive hypervolume)
	// Second item is updated later
	std::vector<std::pair<unsigned int, double> > discarded_islanders;

	std::vector<std::pair<unsigned int, double> > point_pairs;
	// index of currently processed point in the point_pair vector.
	// Initiated to its size (=0) in order to enforce the initial computation on penultimate front.
	unsigned int current_point = point_pairs.size();

	// Stops when we reduce the augmented population to the size of the original population or when the number of discarded islanders reaches the limit
	while (processed_individuals < filtered_immigrants.size() && discarded_islanders.size() < rate_limit) {

		// if current front was exhausted, load next one
		if (current_point == point_pairs.size()) {
			--front_idx;

			// Compute contributions
			std::vector<double> c;

			// If there exist a dominated front for front at index front_idx
			if (front_idx + 1 < fronts_f.size()) {
				std::vector<fitness_vector> merged_front;
				// Reserve the memory and copy the fronts
				merged_front.reserve(fronts_f[front_idx].size() + fronts_f[front_idx + 1].size());

				copy(fronts_f[front_idx].begin(), fronts_f[front_idx].end(), back_inserter(merged_front));
				copy(fronts_f[front_idx + 1].begin(), fronts_f[front_idx +1].end(), back_inserter(merged_front));

				hypervolume hv(merged_front, false);
				c = hv.contributions(refpoint);
			} else {
				hypervolume hv(fronts_f[front_idx], false);
				c = hv.contributions(refpoint);
			}

			// Initiate the pairs and sort by second item (exclusive volume)
			point_pairs.resize(fronts_f[front_idx].size());
			for(unsigned int i = 0 ; i < fronts_f[front_idx].size() ; ++i) {
				point_pairs[i] = std::make_pair(i, c[i]);
			}
			current_point = 0;
			std::sort(point_pairs.begin(), point_pairs.end(), sort_point_pairs_asc);
		}

		unsigned int orig_lc_idx = fronts_i[front_idx][point_pairs[current_point].first];

		if (orig_lc_idx < dest.size()) {
			discarded_islanders.push_back(std::make_pair(orig_lc_idx, 0.0));
		} else {
			++no_discarded_immigrants;
		}

		// Flag given individual as discarded
		g_orig_indices[orig_lc_idx] = 0;

		++processed_individuals;
		++current_point;
	}

	// Number of non-discarded immigrants
	unsigned int no_available_immigrants = boost::numeric_cast<unsigned int>(filtered_immigrants.size() - no_discarded_immigrants);

	// Pairs of (immigrant index, immigrant exclusive hypervolume)
	// Second item is updated later
	std::vector<std::pair<unsigned int, double> > available_immigrants;
	available_immigrants.reserve(no_available_immigrants);
	for(unsigned int idx = dest.size() ; idx < pop_copy.size() ; ++idx) {
		// If the immigrant was not discarded add it to the available set
		if ( g_orig_indices[idx] == 1 ) {
			available_immigrants.push_back(std::make_pair(idx, 0.0));
		}
	}

	// Aggregate all points to establish the hypervolume contribution of available immigrants and discarded islanders
	std::vector<fitness_vector> merged_fronts;
	merged_fronts.reserve(pop_copy.size());

	for(unsigned int idx = 0 ; idx < pop_copy.size() ; ++idx) {
		merged_fronts.push_back(pop_copy.get_individual(idx).cur_f);
	}

	hypervolume hv(merged_fronts, false);
	std::vector<std::pair<unsigned int, double> >::iterator it;

	for(it = available_immigrants.begin() ; it != available_immigrants.end() ; ++it) {
		(*it).second = hv.exclusive((*it).first, refpoint);
	}

	for(it = discarded_islanders.begin() ; it != discarded_islanders.end() ; ++it) {
		(*it).second = hv.exclusive((*it).first, refpoint);
	}

	// Sort islanders and immigrants according to exclusive hypervolume
	sort(available_immigrants.begin(), available_immigrants.end(), hv_fair_r_policy::ind_cmp);
	sort(discarded_islanders.begin(), discarded_islanders.end(), hv_fair_r_policy::ind_cmp);

	// Number of exchanges is the minimum of the number of non discarded immigrants and the number of discarded islanders
	unsigned int no_exchanges = std::min(boost::numeric_cast<unsigned int>(available_immigrants.size()), boost::numeric_cast<unsigned int>(discarded_islanders.size()));

	it = available_immigrants.begin();
	std::vector<std::pair<unsigned int, double> >::reverse_iterator r_it = discarded_islanders.rbegin();

	// Match the best immigrant (forward iterator) with the worst islander (reverse iterator) no_exchanges times.
	for(unsigned int i = 0 ; i < no_exchanges ; ++i) {
		// Break if any islander is better than an immigrant
		if ((*r_it).second > (*it).second) {
			break;
		}
		// Push the pair (islander_idx, fixed_immigrant_idx) to the results
		result.push_back(std::make_pair((*r_it).first, original_immigrant_indices[(*it).first - dest.size()]));
		++r_it;
		++it;
	}

	return result;
}

/// Comparator function used for the sorting of the augmented population fronts
bool hv_fair_r_policy::sort_point_pairs_asc(const std::pair<unsigned int, double> &a, const std::pair<unsigned int, double> &b)
{
	return a.second < b.second;
}

/// Individuals comparator function, used for sorting the individuals according to their exclusive hypervolume (second item)
bool hv_fair_r_policy::ind_cmp(const std::pair<unsigned int, double> &a, const std::pair<unsigned int, double> &b)
{
	return a.second > b.second;
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::migration::hv_fair_r_policy)
