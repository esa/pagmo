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
#include "base.h"
#include <algorithm>

namespace pagmo { namespace util { namespace hv_algorithm {

// comparator function for initial sorting
bool wfg::cmp_points(const fitness_vector &a, const fitness_vector &b) {
	for(int i=a.size() - 1; i >=0 ; --i){
		if (a[i] > b[i]) {
			return true;
		} else if(a[i] < b[i]) {
			return false;
		}
	}
	return false;
};

// Constructor
wfg::wfg(const unsigned int stop_dimension) : m_current_slice(0), m_stop_dimension(stop_dimension) {
	if (stop_dimension < 2 ) {
		pagmo_throw(value_error, "Stop dimension for WFG must be greater than or equal to 2");
	}
}

/// Compute hypervolume 
/**
 * @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume.
 */
double wfg::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	// copy the initial set
	std::vector<fitness_vector> points_cpy(points.begin(), points.end());

	// variable holding the current "depth" of dimension slicing
	// we slice by reducing dimensions from the end
	m_current_slice = r_point.size();

	return compute_hv(points_cpy, r_point);
}

std::vector<fitness_vector> wfg::limitset(const std::vector<fitness_vector> & points, const unsigned int p_idx) const {
	std::vector<fitness_vector> q;
	q.reserve(points.size());

	const fitness_vector& p = points[p_idx];

	for(std::vector<fitness_vector>::size_type idx = p_idx + 1; idx < points.size(); ++idx) {

		fitness_vector s(points[idx]);
		for(fitness_vector::size_type f_idx = 0; f_idx < m_current_slice; ++f_idx) {
			s[f_idx] = fmax(s[f_idx], p[f_idx]);
		}

		int cmp_results[q.size()];

		bool keep_s = true;

		// check whether any point is dominating the point 's'
		for(std::vector<fitness_vector>::size_type q_idx = 0; q_idx < q.size(); ++q_idx) {
			cmp_results[q_idx] = base::dom_cmp(s, q[q_idx], m_current_slice);
			if (cmp_results[q_idx] == 1){
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
				if( cmp_results[cmp_index] == 2 || cmp_results[cmp_index] == 3) {
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

double wfg::compute_hv(const std::vector<fitness_vector> &points, const fitness_vector &r) const {

	std::vector<fitness_vector> points_cpy(points.size());

	for(unsigned int i = 0 ; i < points_cpy.size() ; ++i){
		points_cpy[i] = fitness_vector(points[i].begin(), points[i].begin() + m_current_slice);
	}
	fitness_vector r_cpy(r.begin(), r.begin() + m_current_slice);

	if (m_current_slice == m_stop_dimension) {
		// if already sliced to dimension at which we use another algorithm
		hypervolume hv = hypervolume(points_cpy, false);
		hv.set_copy_points(false);
		return hv.compute(r_cpy);
	} else {
		// else sort the points in preparation for wfg
		sort(points_cpy.begin(), points_cpy.end(), wfg::cmp_points);
	}

	double H = 0.0;
	--m_current_slice;
	for(std::vector<fitness_vector>::size_type i = 0 ; i < points.size() ; ++i) {
		H += fabs((points_cpy[i][m_current_slice] - r_cpy[m_current_slice]) * exclusive_hv(points_cpy, i, r_cpy));
	}
	++m_current_slice;
	return H;
}

double wfg::exclusive_hv(const std::vector<fitness_vector> &points, const unsigned int p_idx, const fitness_vector &r) const {
	std::vector<fitness_vector> q = limitset(points, p_idx);

	double H = base::volume_between(points[p_idx], r, m_current_slice);

	if (q.size() == 1) {
		H -= base::volume_between(q[0], r, m_current_slice );
	} else if (q.size() > 1) {
		H -= compute_hv(q, r);
	}

	return H;
}

// verify_before_compute
/**
 * Verifies whether given algorithm suits the requested data.
 *
 * @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
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

/// Algorithm name
std::string wfg::get_name() const {
	return "WFG algorithm";
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::wfg);
