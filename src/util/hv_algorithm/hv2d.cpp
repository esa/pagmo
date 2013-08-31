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


#include "hv2d.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// Constructor
hv2d::hv2d(const bool initial_sorting) : m_initial_sorting(initial_sorting) { }

/// Computes hypervolume method.
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
double hv2d::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point)
{
	if (points.size() == 0) {
		return 0.0;
	} else if (points.size() == 1) {
		return base::volume_between(points[0], r_point);
	}

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

/// Comparison function for arrays of double.
/**
 * Required by the hv2d::compute method for the sorting of arrays of double*.
 */
bool hv2d::cmp_double_2d(double* a, double* b)
{
	return a[1] < b[1];
}

/// Compute hypervolume method.
/**
 * This method should be used both as a solution to 2-dimensional cases, and as a general termination method for algorithms that reduce n-dimensional problem to 2d.
 * This method is overloaded to work with arrays of double, in order to provide other algorithms that internally work with arrays (such as hv_algorithm::wfg) with an efficient computation.
 *
 * Computational complexity: n*log(n)
 *
 * @param[in] points array of 2-dimensional points
 * @param[in] n_points number of points
 * @param[in] r_point 2-dimensional reference point for the points
 *
 * @return hypervolume
 */
double hv2d::compute(double** points , unsigned int n_points, double* r_point)
{
	if (n_points == 0) {
		return 0.0;
	}
	else if (n_points == 1) {
		return base::volume_between(points[0], r_point, 2);
	}

	if (m_initial_sorting) {
		std::sort(points, points + n_points, hv2d::cmp_double_2d);
	}

	double hypervolume = 0.0;
	for(unsigned int idx = 0; idx < n_points - 1 ; ++idx) {
		double area = (points[idx][1] - points[idx+1][1]) * (points[idx][0] - r_point[0]);
		hypervolume += fabs(area);
	}
	double* last = points[n_points - 1];
	hypervolume += fabs((r_point[1] - last[1]) * (r_point[0] - last[0]));

	return hypervolume;
}

/// Comparison function for sorting of pairs (point, index)
/**
 * Required for hv2d::extreme_contributor method for keeping track of the original indices when sorting.
 */
bool point_pairs_cmp(const std::pair<fitness_vector, unsigned int> &a, const std::pair<fitness_vector, unsigned int> &b)
{
	return a.first[1] > b.first[1];
}

/// Least contributor method
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
unsigned int hv2d::least_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point)
{
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
unsigned int hv2d::greatest_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point)
{
	return extreme_contributor(points, r_point, base::cmp_greatest);
}

/// Extreme contributor method
/**
 * Returns the extreme contributor (least or greatest) for 2d case, depending on the provided comparison function 'cmp_func'
 */
unsigned int hv2d::extreme_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point, bool (*cmp_func)(double, double))
{
	if (points.size() == 1) {
		return 0;
	}
	// introduce a pair of <point, index> in order to return the correct original index
	// otherwise, we would lose the information about the index after the sorting
	std::vector<std::pair<fitness_vector, unsigned int> > points_cpy;
	points_cpy.resize(points.size());
	for (std::vector<std::pair<fitness_vector, unsigned int> >::size_type idx = 0; idx < points.size() ; ++idx) {
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
base_ptr hv2d::clone() const
{
	return base_ptr(new hv2d(*this));
}

/// Algorithm name
std::string hv2d::get_name() const
{
	return "hv2d algorithm";
}

/// Verify input method.
/**
 * Verifies whether the requested data suits the hv2d algorithm.
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the vector of points
 *
 * @throws value_error when trying to compute the hypervolume for the dimension other than 3 or non-maximal reference point
 */
void hv2d::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	if (r_point.size() != 2) {
		pagmo_throw(value_error, "Algorithm hv2d works only for 2-dimensional cases.");
	}

	base::assert_minimisation(points, r_point);
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::hv2d);
