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


#include "wfg.h"

namespace pagmo { namespace util { namespace hv_algorithm {


/// Compute hypervolume 
/**
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume.
 */
double wfg::compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	// copy the initial set
	std::vector<fitness_vector> points_cpy(points.begin(), points.end());

	// sort the initial set by first dimension
	sort(points_cpy.begin(), points_cpy.end(), fitness_vector_cmp(0,'<'));

	return compute_hv(points_cpy, r_point);
}

std::vector<fitness_vector> wfg::limitset(const std::vector<fitness_vector> & points, unsigned int p_idx) {
	std::vector<fitness_vector> q;

	for(std::vector<fitness_vector>::size_type idx = p_idx + 1; idx < points.size(); ++idx) {
		fitness_vector s(points[idx]);
		for(fitness_vector::size_type f_idx = 0; f_idx < points[idx].size(); ++f_idx) {
			s[f_idx] = fmax(s[f_idx], points[p_idx][f_idx]);
		}

		int cmp_results[q.size()];

		bool keep_s = true;

		// check whether any point is dominating the point 's'
		for(std::vector<fitness_vector>::size_type q_idx = 0; q_idx < q.size(); ++q_idx) {
			cmp_results[q_idx] = dom_cmp(s,q[q_idx]);
			if (cmp_results[q_idx] == -1){
				keep_s = false;
				break;
			}
		}
		// if neither is, remove points dominated by 's' (we store that during the first loop)
		if( keep_s ) {
			std::vector<fitness_vector>::iterator q_it = q.begin();
			int cmp_index = 0;
			// loop over points by iterator, but maintain the counter of previously stored comparison results
			while(q_it != q.end()) {
				if( cmp_results[cmp_index] == 1) {
					q_it = q.erase(q_it);
				} else {
					++q_it;
				}
				++cmp_index;
			}
			q.push_back(s);
		}
	}
	return q;
}

double wfg::inclusive_hv(const fitness_vector &p, const fitness_vector &r) {
	double total_hv = 1.0;
	for(fitness_vector::size_type idx = 0 ; idx < p.size() ; ++idx) {
		total_hv *= (r[idx] - p[idx]);
	}
	return fabs(total_hv);
}


double wfg::compute_hv(const std::vector<fitness_vector> &points, const fitness_vector &r) {
	double H = 0.0;
	for(std::vector<fitness_vector>::size_type idx = 0 ; idx < points.size() ; ++idx) {
		H += exclusive_hv(points, idx, r);
	}
	return H;
}

double wfg::exclusive_hv(const std::vector<fitness_vector> &points, unsigned int p_idx, const fitness_vector &r) {
	std::vector<fitness_vector> q = limitset(points, p_idx);
	return inclusive_hv(points[p_idx], r) - compute_hv(q, r);
}

/// dominance comparison
/** 
 * Establishes the relationship between two points.
 *
 * return -1 if 'a' IS DOMINATED BY 'b'
 * return 1 if ('a' DOMINATES 'b') or ('a' EQUAL TO 'b')
 * return 0 otherwise
 */
int wfg::dom_cmp(const fitness_vector &a, const fitness_vector &b) {
	for(fitness_vector::size_type i = 0; i < a.size() ; ++i) {
		if (a[i] > b[i]) {
			for(fitness_vector::size_type j = i + 1; j < a.size() ; ++j) {
				if (a[j] < b[j]) {
					return 0;
				}
			}
			return -1;
		}
		else if (a[i] > b[i]) {
			for(fitness_vector::size_type j = i + 1 ; j < a.size() ; ++j) {
				if (a[j] < b[j]) {
					return 0;
				}
			}
			return 1;
		}
	}
	return 1;
}

// verify_before_compute
/**
 * Verifies whether given algorithm suits the requested data.
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the vector of points
 *
 * @throws value_error when trying to compute the hypervolume for the non-maximal reference point
 */
void wfg::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	base::assert_maximal_reference_point(points, r_point);
}

/// Clone method.
base_ptr wfg::clone() const
{
	return base_ptr(new wfg(*this));
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::wfg);
