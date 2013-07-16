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


#include "native2d.h"

namespace pagmo { namespace util { namespace hv_algorithm {

// Computes hypervolume method.
/**
 * This method should be used both as a solution to 2D cases, and as a general termination method for algorithms that reduce n-dimensional problem to 2D.
 *
 * Computational complexity: n*log(n)
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume
 */
double native2d::compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point)
{
	std::vector<fitness_vector> points_cpy(points.begin(), points.end());
	sort(points_cpy.begin(), points_cpy.end(), fitness_vector_cmp(0, '<'));
	double hypervolume = 0.0;
	for(std::vector<fitness_vector>::size_type idx = 0; idx < points_cpy.size() - 1 ; ++idx) {
		double area = (points_cpy[idx][0] - points_cpy[idx+1][0]) * (points_cpy[idx][1] - r_point[1]);
		hypervolume += fabs(area);
	}
	fitness_vector &last = points_cpy.back();
	hypervolume += fabs((r_point[0] - last[0]) * (r_point[1] - last[1]));

	return hypervolume;
}

// custom comparison method for sorting pairs of (point, index)
// required for native2d::least_contributor method
bool point_pairs_cmp(const std::pair<fitness_vector, unsigned int> &a, const std::pair<fitness_vector, unsigned int> &b) {
	return a.first[0] > b.first[0];
}

// Least contributing point method
/**
 * This method overloads the default method for calculating the least contributing point.
 * It uses a similar method as compute method, that is, sort each point by one dimension, and then find the least contributor in a linear fashion.
 *
 * @see Nicola Beume, Boris Naujoks, Michael Emmerich, "SMS-EMOA: Multiobjective selection based on dominated hypervolume", Section 2.3. "Calculation of contributing hypervolume".
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return index of the least contributing point
 */
unsigned int native2d::least_contributor(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	if (points.size() == 1) {
		return 0;
	}
	// introduce a pair of <point, index> in order to return the correct original index
	// otherwise, we would loose index information after sorting
	std::vector<std::pair<fitness_vector, unsigned int> > points_cpy;
	points_cpy.resize(points.size());
	for(std::vector<std::pair<fitness_vector, unsigned int> >::size_type idx = 0; idx < points.size() ; ++idx) {
		points_cpy[idx] = std::pair<fitness_vector, unsigned int>(points[idx], idx);
	}

	sort(points_cpy.begin(), points_cpy.end(), point_pairs_cmp);

	// compute first point separately
	double least_contrib_idx = 0;
	double least_hv = fabs((points_cpy[0].first[0] - r_point[0]) * (points_cpy[0].first[1] - points_cpy[1].first[1]));

	// compute points 2nd to (m-1)th
	for(std::vector<fitness_vector>::size_type idx = 1; idx < points_cpy.size() - 1 ; ++idx) {
		double exclusive_hv = fabs((points_cpy[idx].first[0] - points_cpy[idx - 1].first[0]) * (points_cpy[idx].first[1] - points_cpy[idx+1].first[1]));
		if (exclusive_hv < least_hv) {
			least_hv = exclusive_hv;
			least_contrib_idx = idx;
		}
	}
	unsigned int last_idx = points_cpy.size() - 1;

	// compute last point separately
	double last_exclusive_hv = fabs((points_cpy[last_idx].first[0] - points_cpy[last_idx - 1].first[0]) * (points_cpy[last_idx].first[1] - r_point[1]));
	if (last_exclusive_hv < least_hv) {
		least_contrib_idx = last_idx;
	}

	return points_cpy[least_contrib_idx].second;
}

/// Clone method.
base_ptr native2d::clone() const
{
	return base_ptr(new native2d(*this));
}

/// Algorithm name
std::string native2d::get_name() const {
	return "Native2D algorithm";
}

/// Verify input method.
/**
 * Verifies whether Native2D algorithm suits the requested data.
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the vector of points
 *
 * @throws value_error when trying to compute the hypervolume for the dimension other than 3 or non-maximal reference point
 */
void native2d::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	if (r_point.size() != 2) {
		pagmo_throw(value_error, "native2d method method works only for 2-dimensional cases.");
	}

	base::assert_maximal_reference_point(points, r_point);
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::native2d);
