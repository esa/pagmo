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

#include "base.h"


namespace pagmo { namespace util { namespace hv_algorithm {

/// Destructor required for pure virtual methods
base::~base() {}

/// Assert that reference point dominates every other point from the set.
/**
 * This is a method that can be referenced from verify_before_compute method.
 * The method checks whether provided reference point is "no worse" and in at least one objective "better" for each of the points from the set.
 *
 * @param[in] points - vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point - distinguished "reference point".
*/
void base::assert_maximal_reference_point(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	for(std::vector<fitness_vector>::size_type idx = 0 ; idx < points.size() ; ++idx) {
		bool outside_bounds = false;
		bool all_equal = true;

		for(fitness_vector::size_type f_idx = 0 ; f_idx < points[idx].size() ; ++f_idx) {
			outside_bounds |= (r_point[f_idx] < points[idx][f_idx]);
			all_equal &= (r_point[f_idx] == points[idx][f_idx]);
		}
		if (all_equal || outside_bounds) {
			std::stringstream ss;
			std::string str_p("("), str_r("(");
			for(fitness_vector::size_type f_idx = 0 ; f_idx < points[idx].size() ; ++f_idx) {
					str_p += boost::lexical_cast<std::string>(points[idx][f_idx]);
					str_r += boost::lexical_cast<std::string>(r_point[f_idx]);
				if (f_idx < points[idx].size() - 1) {
					str_p += ", ";
					str_r += ", ";
				} else {
					str_p += ")";
					str_r += ")";
				}
			}
			ss << "Reference point is invalid: another point seems to be outside the reference point boundary, or be equal to it:" << std::endl;
			ss << " P[" << idx << "]\t= " << str_p << std::endl;
			ss << " R\t= " << str_r << std::endl;
			pagmo_throw(value_error, ss.str());
		}
	}
}

// Compute volume between two points
/**
 * Calculates the volume between points a and b (as defined for n-dimensional Euclidean spaces).
 *
 * @param[in] a first point defining the hypercube
 * @param[in] b second point defining the hypercube
 *
 * @return volume of hypercube defined by points a and b
 */
double base::volume_between(const fitness_vector &a, const fitness_vector &b) const {
	pagmo_assert(a.size() == b.size());
	pagmo_assert(a.size() >= 2);
	double volume = 1.0;
	for (fitness_vector::size_type idx = 0; idx < a.size() ; ++idx) {
		volume *= (a[idx] - b[idx]);
	}
	return (volume < 0 ? -volume : volume);
}

/// Exclusive hypervolume method
/**
 * This method computes the exclusive hypervolume for given individual.
 * It accepts a list of points as an input, and the distinguished "reference point".
 * Hypervolume is then computed as a joint hypervolume of hypercubes, generated pairwise with the reference point and each point from the set.
 *
 * @param[in] p_idx index of the individual
 * @param[in] points vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point distinguished "reference point".
 *
 * @return exlusive hypervolume contributed by the individual at index p_idx
 */
double base::exclusive(const unsigned int p_idx, const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	double hypvol_total = compute(points, r_point);
	if (points.size() == 1) {
		return hypvol_total;
	}
	std::vector<fitness_vector> points_less;
	points_less.reserve(points.size()-1);
	copy(points.begin(), points.begin()+p_idx, back_inserter(points_less));
	copy(points.begin()+p_idx+1, points.end(), back_inserter(points_less));
	return hypvol_total - compute(points_less, r_point);
}

/// Least contributing point method
/**
 * This method establishes the individual that contributes the least to the hypervolume.
 * By default it computes each individual contribution, and chooses the one that contributes the least.
 * Other algorithms may overload this method for a more efficient solution.
 *
 * @param[in] points vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point distinguished "reference point".
 *
 * @return index of the least contributing point
 */
unsigned int base::least_contributor(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	if (points.size() == 1) {
		return 0;
	}

	double hv_total = compute(points, r_point);
	std::vector<fitness_vector> points_cpy(points.begin() + 1, points.end());
	double idx_min = 0;
	double hv_min = hv_total - compute(points_cpy, r_point);

	for(unsigned int idx = 1 ; idx < points.size() ; ++idx) {
		std::vector<fitness_vector> points_less;
		points_less.reserve(points.size()-1);
		copy(points.begin(), points.begin()+idx, back_inserter(points_less));
		copy(points.begin()+idx+1, points.end(), back_inserter(points_less));
		double hv = hv_total - compute(points_less, r_point);
		if (hv < hv_min) {
			hv_min = hv;
			idx_min = idx;
		}
	}

	return idx_min;
}


///Constructor of the comparator object
/**
 * Create a comparator object, that compares items by given dimension, according to given inequality function.
 *
 * @param[in] dim dimension index by which we compare the fitness vectors
 * @param[in] cmp_type inequality expression used for comparison, either character '<' or '>'
 */
fitness_vector_cmp::fitness_vector_cmp(int dim, char cmp_type) {
	if (cmp_type == '<') {
		m_cmp_obj = boost::shared_ptr<cmp_fun>( new cmp_le(dim));
	}
	else {
		m_cmp_obj = boost::shared_ptr<cmp_fun>( new cmp_ge(dim));
	}
}

/// Get algorithm's name.
/**
 * Default implementation will return the algorithm's mangled C++ name.
 *
 * @return name of the algorithm.
 */
std::string base::get_name() const
{
	return typeid(*this).name();
}

} } }
