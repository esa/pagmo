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

#include "base.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// Destructor required for pure virtual methods
base::~base() {}

/// Assert that reference point dominates every other point from the set.
/**
 * This is a method that can be referenced from verify_before_compute method.
 * The method checks whether the provided reference point fits the minimisation assumption, e.g.,
 * reference point must be "no worse" and in at least one objective and "better" for each of the points from the set.
 *
 * @param[in] points - vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point - distinguished "reference point".
*/
void base::assert_minimisation(const std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	for(std::vector<fitness_vector>::size_type idx = 0 ; idx < points.size() ; ++idx) {
		bool outside_bounds = false;
		bool all_equal = true;

		for(fitness_vector::size_type f_idx = 0 ; f_idx < points[idx].size() ; ++f_idx) {
			outside_bounds |= (r_point[f_idx] < points[idx][f_idx]);
			all_equal &= (r_point[f_idx] == points[idx][f_idx]);
		}
		if (all_equal || outside_bounds) {
			// Prepare error message.
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
double base::exclusive(const unsigned int p_idx, std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	if (points.size() == 1) {
		return compute(points, r_point);
	}
	std::vector<fitness_vector> points_less;
	points_less.reserve(points.size() - 1);
	copy(points.begin(), points.begin() + p_idx, back_inserter(points_less));
	copy(points.begin() + p_idx + 1, points.end(), back_inserter(points_less));

	return compute(points, r_point) - compute(points_less, r_point);
}

/// compute the extreme contributor
/**
 * Computes the index of the individual that contributes the most or the least to the hypervolume (depending on the prodivded comparison function)
 */
unsigned int base::extreme_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point, bool (*cmp_func)(double, double)) const
{
	// Trivial case
	if (points.size() == 1) {
		return 0;
	}

	std::vector<double> c = contributions(points, r_point);

	unsigned int idx_extreme = 0;

	// Check the remaining ones using the provided comparison function
	for(unsigned int idx = 1 ; idx < c.size() ; ++idx) {
		if (cmp_func(c[idx], c[idx_extreme])) {
			idx_extreme = idx;
		}
	}

	return idx_extreme;
}

/// Comparison method for the least contributor
/**
 * This method is used as a comparison function for the extreme contributor method which may be overloaded by hv algorithms.
 * In such case, this method can determine, given two contributions of points, which one is the "smaller" contributor.
 *
 * @param[in] a first contribution of a point
 * @param[in] b second contribution of a point
 *
 * @return true if contribution 'a' is lesser than contribution 'b'
 */
bool base::cmp_least(const double a, const double b)
{
	return a < b;
}


/// Comparison method for the least contributor
/**
 * This method is used as a comparison function for the extreme contributor method which may be overloaded by hv algorithms.
 * In such case, this method can determine, given two contributions of points, which one is the "greater" contributor.
 *
 * @param[in] a first contribution of a point
 * @param[in] b second contribution of a point
 *
 * @return true if contribution 'a' is greater than contribution 'b'
 */
bool base::cmp_greatest(const double a, const double b)
{
	return a > b;
}

/// Least contributor method
/**
 * This method establishes the individual that contributes the least to the hypervolume.
 * By default it computes each individual contribution, and chooses the one with the lowest contribution.
 * Other algorithms may overload this method for a more efficient solution.
 *
 * @param[in] points vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point distinguished "reference point".
 *
 * @return index of the least contributor
 */
unsigned int base::least_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	return extreme_contributor(points, r_point, base::cmp_least);
}

/// Greatest contributor method
/**
 * This method establishes the individual that contributes the most to the hypervolume.
 * By default it computes each individual contribution, and chooses the one with the highest contribution.
 * Other algorithms may overload this method for a more efficient solution.
 *
 * @param[in] points vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point distinguished "reference point".
 *
 * @return index of the greatest contributor
 */
unsigned int base::greatest_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	return extreme_contributor(points, r_point, base::cmp_greatest);
}

/// Contributions method
/**
 * This methods return the exclusive contribution to the hypervolume for each point.
 * Main reason for this method is the fact that in most cases the explicit request for all contributions
 * can be done more efficiently (may vary depending on the provided hv_algorithm) than executing "exclusive" method in a loop.
 *
 * This base method uses a very naive approach, which in fact is only slightly more efficient than calling "exclusive" method successively.
 *
 * @param[in] points vector of fitness_vectors for which the contributions are computed
 * @param[in] r_point distinguished "reference point".
 * @return vector of exclusive contributions by every point
 */
std::vector<double> base::contributions(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	std::vector<double> c;
	c.reserve(points.size());

	// Trivial case
	if (points.size() == 1) {
		c.push_back(volume_between(points[0], r_point));
		return c;
	}

	// Compute the total hypervolume for the reference
	std::vector<fitness_vector> points_cpy(points.begin(), points.end());
	double hv_total = compute(points_cpy, r_point);

	// Points[0] as a first candidate
	points_cpy = std::vector<fitness_vector>(points.begin() + 1, points.end());
	c.push_back(hv_total - compute(points_cpy, r_point));

	// Check the remaining ones using the provided comparison function
	for(unsigned int idx = 1 ; idx < points.size() ; ++idx) {
		std::vector<fitness_vector> points_less;
		points_less.reserve(points.size() - 1);
		copy(points.begin(), points.begin() + idx, back_inserter(points_less));
		copy(points.begin() + idx + 1, points.end(), back_inserter(points_less));
		double delta = hv_total - compute(points_less, r_point);

		if (fabs(delta) < 1e-8) {
			delta = 0.0;
		}
		c.push_back(delta);
	}

	return c;
}

/// Compute volume between two points
/**
 * Calculates the volume between points a and b (as defined for n-dimensional Euclidean spaces).
 *
 * @param[in] a first point defining the hypercube
 * @param[in] b second point defining the hypercube
 * @param[in] dim_bound dimension boundary for the volume. If equal to 0 (default value), then compute the volume of whole vector. Any positive number limits the computation from dimension 1 to dim_bound INCLUSIVE.
 *
 * @return volume of hypercube defined by points a and b
 */
double base::volume_between(const fitness_vector &a, const fitness_vector &b, unsigned int dim_bound)
{
	if (dim_bound == 0) {
		dim_bound = a.size();
	}
	double volume = 1.0;
	for (fitness_vector::size_type idx = 0; idx < dim_bound ; ++idx) {
		volume *= (a[idx] - b[idx]);
	}
	return (volume < 0 ? -volume : volume);
}

/// Compute volume between two points
/**
 * Calculates the volume between points a and b (as defined for n-dimensional Euclidean spaces).
 *
 * @param[in] a first point defining the hypercube
 * @param[in] b second point defining the hypercube
 * @param[in] size dimension of the vectors.
 *
 * @return volume of hypercube defined by points a and b
 */
double base::volume_between(double* a, double* b, unsigned int size)
{
	double volume = 1.0;
	while(size--) {
		volume *= (b[size] - a[size]);
	}
	return (volume < 0 ? -volume : volume);
}

/// Dominance comparison method
/**
 * Establishes the domination relationship between two points.
 *
 * returns DOM_CMP_B_DOMINATES_A if point 'b' DOMINATES point 'a'
 * returns DOM_CMP_A_DOMINATES_B if point 'a' DOMINATES point 'b'
 * returns DOM_CMP_A_B_EQUAL if point 'a' IS EQUAL TO 'b'
 * returns DOM_CMP_INCOMPARABLE otherwise
 */
int base::dom_cmp(const fitness_vector &a, const fitness_vector &b, unsigned int dim_bound)
{
	if (dim_bound == 0) {
		dim_bound = a.size();
	}
	for(fitness_vector::size_type i = 0; i < dim_bound ; ++i) {
		if (a[i] > b[i]) {
			for(fitness_vector::size_type j = i + 1; j < dim_bound ; ++j) {
				if (a[j] < b[j]) {
					return DOM_CMP_INCOMPARABLE;
				}
			}
			return DOM_CMP_B_DOMINATES_A;
		}
		else if (a[i] < b[i]) {
			for(fitness_vector::size_type j = i + 1 ; j < dim_bound ; ++j) {
				if (a[j] > b[j]) {
					return DOM_CMP_INCOMPARABLE;
				}
			}
			return DOM_CMP_A_DOMINATES_B;
		}
	}
	return DOM_CMP_A_B_EQUAL;
}

/// Dominance comparison method
/**
 * Establishes the domination relationship between two points (overloaded for double*);
 *
 * returns DOM_CMP_B_DOMINATES_A if point 'b' DOMINATES point 'a'
 * returns DOM_CMP_A_DOMINATES_B if point 'a' DOMINATES point 'b'
 * returns DOM_CMP_A_B_EQUAL if point 'a' IS EQUAL TO 'b'
 * returns DOM_CMP_INCOMPARABLE otherwise
 */
int base::dom_cmp(double* a, double* b, unsigned int size)
{
	for(fitness_vector::size_type i = 0; i < size ; ++i) {
		if (a[i] > b[i]) {
			for(fitness_vector::size_type j = i + 1; j < size ; ++j) {
				if (a[j] < b[j]) {
					return DOM_CMP_INCOMPARABLE;
				}
			}
			return DOM_CMP_B_DOMINATES_A;
		}
		else if (a[i] < b[i]) {
			for(fitness_vector::size_type j = i + 1 ; j < size ; ++j) {
				if (a[j] > b[j]) {
					return DOM_CMP_INCOMPARABLE;
				}
			}
			return DOM_CMP_A_DOMINATES_B;
		}
	}
	return DOM_CMP_A_B_EQUAL;
}

///Constructor of the comparator object
/**
 * Create a comparator object, that compares items by given dimension, according to given inequality function.
 *
 * @param[in] dim dimension index by which we compare the fitness vectors
 * @param[in] cmp_type inequality expression used for comparison, either character '<' or '>'
 */
fitness_vector_cmp::fitness_vector_cmp(int dim, char cmp_type)
{
	if (cmp_type == '<') {
		m_cmp_obj = boost::shared_ptr<cmp_fun>(new cmp_le(dim));
	}
	else {
		m_cmp_obj = boost::shared_ptr<cmp_fun>(new cmp_ge(dim));
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
