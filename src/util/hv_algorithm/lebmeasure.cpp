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

#include "lebmeasure.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// lebmeasure::compute
/**
 * Computes hypervolume indicator for given pareto set using the LebMeasure algorithm.
 *
 * @param[in] points set of points describing the hypervolume
 * @param[in] r_point reference point for the hypervolume
 *
 * @throws value_error if sizes of reference point and pareto set do not match, when the pareto set is empty and when the dimension is lesser than 3
 *
 * @return hypervolume of the pareto set.
 */
double lebmeasure::compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	lebmeasure_points point_set;
	fitness_vector::size_type f_dim = points[0].size();
	for (std::vector<fitness_vector>::size_type idx = 0 ; idx < points.size() ; ++idx) {
		point_set.push_back(std::make_pair(points[idx], f_dim));
	}

	double hypervolume = 0.0;
	while(point_set.size() > 0) {
		std::pair<fitness_vector, fitness_vector::size_type> head = point_set.front();
		fitness_vector p = head.first;
		fitness_vector::size_type z = head.second;
		point_set.pop_front();
		fitness_vector a = get_opposite_point(p, point_set, r_point);
		hypervolume += volume_between(p, a);
		lebmeasure_points ql = generate_spawns(p, z, a, point_set, r_point);
		point_set.insert(point_set.end(), ql.begin(), ql.end());
	}
	return hypervolume;
}

// lebmeasure::verify_before_compute
/**
 * Verifies whether this method is aplicable to given hypervolume problem.
 *
 * @param[in] points set of points describing the hypervolume
 * @param[in] r_point reference point for the hypervolume
 *
 * @throws value_error if the reference point is of the the dimension lesser than 3 or non-maximal reference point
 */
void lebmeasure::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	if (r_point.size() < 3) {
		pagmo_throw(value_error, "Hypervolume of dimension lesser than 3 is not allowed for this method, use optimal 2D instead");
	}
	base::assert_maximal_reference_point(points, r_point);
}


/// Clone method.
base_ptr lebmeasure::clone() const
{
	return base_ptr(new lebmeasure(*this));
}

/// Dominated method
/**
 * Determinines whether given point p is strictly dominated by any point from given set of points.
 *
 * @param[in] p fitness vector describing given point in space
 * @param[in] point_set set of pairs <point, spawn_dimension>, p is compared against for strict dominance
 *
 * @return true when p is strictly dominated by any point from point_set, false otherwise.
 */
bool lebmeasure::dominated(const fitness_vector &p, lebmeasure_points &point_set)
{
	for (lebmeasure_points::iterator it = point_set.begin() ; it != point_set.end() ; ++it) {
		bool is_dominated = true;
		for (fitness_vector::size_type idx_dim = 0 ; idx_dim < p.size() ; ++idx_dim) {
			if ((*it).first[idx_dim] > p[idx_dim]) {
				is_dominated = false;
				break;
			}
		}
		if (is_dominated) {
			return true;
		}
	}
	return false;
}

// Generates "spawn" points for given set of points.
/**
 * @param[in] p fitness vector describing given point in fitness space, from which spawn points are generated
 * @param[in] z number of dimensions to consider
 * @param[in] a point opposite to p (see method "get_opposite_point")
 * @param[in] point_set set of pairs (point, spawn_dimension) forcing some restrictions on final spawned point set
 * @param[in] r_point reference point for hypervolume
 *
 * @return vector of points spawned from p along with the dimensions that spawned given point.
 */
lebmeasure_points lebmeasure::generate_spawns(const fitness_vector &p, const fitness_vector::size_type z, const fitness_vector &a, lebmeasure_points &point_set, const fitness_vector &r_point)
{
	lebmeasure_points ql; //pair here
	for (fitness_vector::size_type idx_dim = 0 ; idx_dim < z ; ++idx_dim) {
		if (a[idx_dim] != r_point[idx_dim]) {
			fitness_vector q = fitness_vector(p.begin(), p.end());
			q[idx_dim] = a[idx_dim];
			if (!dominated(q, point_set))
				ql.push_back(std::make_pair(q, idx_dim+1));
		}
	}
	return ql;
}

// Compute volume between two points
/**
 * Calculates the volume between points a and b (as defined for n-dimensional Euclidian spaces).
 *
 * @param[in] a first point defining the hypercube
 * @param[in] b second point defining the hypercube
 *
 * @return volume of hypercube defined by points a and b
 */
double lebmeasure::volume_between(const fitness_vector &a, const fitness_vector &b)
{
	pagmo_assert(a.size() == b.size());
	pagmo_assert(a.size() >= 2);
	double volume = 1.0;
	for (fitness_vector::size_type idx = 0; idx < a.size() ; ++idx) {
		volume *= (a[idx] - b[idx]);
	}
	return (volume < 0 ? -volume : volume);
}

// Get point opposite to p
/**
 * Generates the point opposite to p, restricted to points from the point_set.
 * It tries to find the point that is closest to original point p, restricted (by each dimension) by the points in the point_set.
 *
 * @param[in] p point for which the opposite point is generated
 * @param[in] point_set set of points forcing restrictions on how "far" the opposite point can be
 * @param[in] r_point reference point of the hypervolume
 *
 * @return point opposite to point p
 */
fitness_vector lebmeasure::get_opposite_point(const fitness_vector &p, lebmeasure_points &point_set, const fitness_vector &r_point)
{
	fitness_vector a(p.size());
	for (fitness_vector::size_type idx_dim = 0 ; idx_dim < p.size() ; ++idx_dim) {
		double next_closer = r_point[idx_dim];
		for(lebmeasure_points::iterator it = point_set.begin() ; it != point_set.end() ; ++it) {
			if ( (*it).first[idx_dim] > p[idx_dim] && (*it).first[idx_dim] < next_closer) {
				next_closer = (*it).first[idx_dim];
			}
		}
		a[idx_dim] = next_closer;
	}
	return a;
}


} } }
