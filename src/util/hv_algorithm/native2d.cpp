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

// Constructor
native2d::native2d(const bool initial_sorting) : m_initial_sorting(initial_sorting) { }

// Computes hypervolume method.
/**
 * This method should be used both as a solution to 2-dimensional cases, and as a general termination method for algorithms that reduce n-dimensional problem to 2D.
 *
 * Computational complexity: n*log(n)
 *
 * @param[in] points vector of points containing the 2-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume
 */
double native2d::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point)
{
	if (points.size() == 0) {
		return 0.0;
	} else if (points.size() == 1) {
		return base::volume_between(points[0], r_point);
	}
	//std::vector<fitness_vector> points_cpy(points.begin(), points.end());

	if (m_initial_sorting) {
		sort(points.begin(), points.end(), fitness_vector_cmp(1, '<'));
	}

	double hypervolume = 0.0;
	for(std::vector<fitness_vector>::size_type idx = 0; idx < points.size() - 1 ; ++idx) {
		double area = (points[idx][1] - points[idx+1][1]) * (points[idx][0] - r_point[0]);
		hypervolume += fabs(area);
	}
	fitness_vector &last = points.back();
	hypervolume += fabs((r_point[1] - last[1]) * (r_point[0] - last[0]));

	return hypervolume;
}

bool native2d::cmp_double_2d(double* a, double* b) {
	return a[1] < b[1];
}

double native2d::compute(double**points , unsigned int n_points, double* r_point)
{
	if (n_points == 0) {
		return 0.0;
	}
	else if (n_points == 1) {
		return base::volume_between(points[0], r_point, 2);
	}

	if (m_initial_sorting) {
		std::sort(points, points + n_points, native2d::cmp_double_2d);
	}

	double hypervolume = 0.0;
	for(unsigned int idx = 0; idx < n_points - 1 ; ++idx) {
		double area = (points[idx][1] - points[idx+1][1]) * (points[idx][0] - r_point[0]);
		hypervolume += fabs(area);
	}
	//fitness_vector &last = points.back();
	double* last = points[n_points - 1];
	hypervolume += fabs((r_point[1] - last[1]) * (r_point[0] - last[0]));

	return hypervolume;
}

// custom comparison method for sorting pairs of (point, index)
// required for native2d::least_contributor method
bool point_pairs_cmp(const std::pair<fitness_vector, unsigned int> &a, const std::pair<fitness_vector, unsigned int> &b) {
	return a.first[1] > b.first[1];
}

// Least contributor method
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
unsigned int native2d::least_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	return extreme_contributor(points, r_point, base::cmp_least);
}

// Greatest contributor method
/**
 * This method overloads the default method for calculating the greatest contributor among the points
 * It uses a similar method as compute method, that is, sort each point by one dimension, and then find the greatest contributor in a linear fashion.
 *
 * @see Nicola Beume, Boris Naujoks, Michael Emmerich, "SMS-EMOA: Multiobjective selection based on dominated hypervolume", Section 2.3. "Calculation of contributing hypervolume".
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return index of the greatest contributor
 */
unsigned int native2d::greatest_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	return extreme_contributor(points, r_point, base::cmp_greatest);
}

// extreme contributor method
// returns the extreme contributor (least or greatest) for 2d case, dependent on the provided comparison function
unsigned int native2d::extreme_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point, bool (*cmp_func)(double, double)) {
	if (points.size() == 1) {
		return 0;
	}
	// introduce a pair of <point, index> in order to return the correct original index
	// otherwise, we would lose the information about the index after the sorting
	std::vector<std::pair<fitness_vector, unsigned int> > points_cpy;
	points_cpy.resize(points.size());
	for(std::vector<std::pair<fitness_vector, unsigned int> >::size_type idx = 0; idx < points.size() ; ++idx) {
		points_cpy[idx] = std::pair<fitness_vector, unsigned int>(points[idx], idx);
	}

	sort(points_cpy.begin(), points_cpy.end(), point_pairs_cmp);

	// compute first point separately
	double extreme_contrib_idx = 0;
	double extreme_hv = fabs((points_cpy[0].first[1] - r_point[1]) * (points_cpy[0].first[0] - points_cpy[1].first[0]));

	// compute points 2nd to (m-1)th
	for(std::vector<fitness_vector>::size_type idx = 1; idx < points_cpy.size() - 1 ; ++idx) {
		double exclusive_hv = fabs((points_cpy[idx].first[1] - points_cpy[idx - 1].first[1]) * (points_cpy[idx].first[0] - points_cpy[idx+1].first[0]));
		if (cmp_func(exclusive_hv, extreme_hv)) {
			extreme_hv = exclusive_hv;
			extreme_contrib_idx = idx;
		}
	}
	unsigned int last_idx = points_cpy.size() - 1;

	// compute last point separately
	double last_exclusive_hv = fabs((points_cpy[last_idx].first[1] - points_cpy[last_idx - 1].first[1]) * (points_cpy[last_idx].first[0] - r_point[0]));
	if (cmp_func(last_exclusive_hv, extreme_hv)) {
		extreme_contrib_idx = last_idx;
	}

	return points_cpy[extreme_contrib_idx].second;
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
		pagmo_throw(value_error, "Algorithm Native2D method works only for 2-dimensional cases.");
	}

	base::assert_maximal_reference_point(points, r_point);
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::native2d);
