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


#include "hv2d.h"
#include "hv3d.h"

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
double hv2d::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
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

	// width of the sweeping line
	double w = r_point[0] - points[0][0];
	for(unsigned int idx = 0 ; idx < points.size() - 1 ; ++idx) {
		hypervolume += (points[idx + 1][1] - points[idx][1]) * w;
		w = std::max(w, r_point[0] - points[idx+1][0]);
	}
	hypervolume += (r_point[1] - points[points.size() - 1][1]) * w;

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
double hv2d::compute(double** points, unsigned int n_points, double* r_point) const
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

	// width of the sweeping line
	double w = r_point[0] - points[0][0];
	for(unsigned int idx = 0 ; idx < n_points - 1 ; ++idx) {
		hypervolume += (points[idx + 1][1] - points[idx][1]) * w;
		w = std::max(w, r_point[0] - points[idx+1][0]);
	}
	hypervolume += (r_point[1] - points[n_points - 1][1]) * w;

	return hypervolume;
}

/// Contributions method
/**
 * Computes the contributions of each point by invoking the HV3D algorithm with mock third dimension.
 *
 * @param[in] points vector of points containing the 2-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 * @return vector of exclusive contributions by every point
 */
std::vector<double> hv2d::contributions(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	std::vector<fitness_vector> new_points(points.size(), fitness_vector(3, 0.0));
	fitness_vector new_r(r_point);
	new_r.push_back(1.0);

	for(unsigned int i = 0 ; i < points.size() ; ++i) {
		new_points[i][0] = points[i][0];
		new_points[i][1] = points[i][1];
		new_points[i][2] = 0.0;
	}
	// Set sorting to off since contributions are sorted by third dimension
	return hv3d(false).contributions(new_points, new_r);
}


/// Comparison function for sorting of pairs (point, index)
/**
 * Required for hv2d::extreme_contributor method for keeping track of the original indices when sorting.
 */
bool hv2d::point_pairs_cmp(const std::pair<fitness_vector, unsigned int> &a, const std::pair<fitness_vector, unsigned int> &b)
{
	return a.first[1] > b.first[1];
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

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::hv2d)
