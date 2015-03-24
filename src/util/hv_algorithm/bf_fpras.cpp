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


#include "bf_fpras.h"
#include <algorithm>

namespace pagmo { namespace util { namespace hv_algorithm {

/// Constructor
/**
 * Constructs an instance of the algorithm
 *
 * @param[in] eps accuracy of the approximation
 * @param[in] delta confidence of the approximation
 */
bf_fpras::bf_fpras(const double eps, const double delta) : m_eps(eps), m_delta(delta) { }

/// Verify before compute
/**
 * Verifies whether given algorithm suits the requested data.
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the vector of points
 *
 * @throws value_error when trying to compute the hypervolume for the non-maximal reference point
 */
void bf_fpras::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	base::assert_minimisation(points, r_point);
}

/// Compute method
/**
 * Compute the hypervolume using FPRAS.
 *
 * @see "Approximating the volume of unions and intersections of high-dimensional geometric objects", Karl Bringmann, Tobias Friedrich.
 *
 * @param[in] points vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point distinguished "reference point".
 *
 * @return approximated hypervolume
 */
double bf_fpras::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	unsigned int n = points.size();
	unsigned int dim = r_point.size();
	boost::uint_fast64_t T = static_cast<boost::uint_fast64_t>( 12. * std::log( 1. / m_delta ) / std::log( 2. ) * n / m_eps / m_eps );

	// Partial sums of consecutive boxes
	std::vector<double> sums(n, 0.0);

	// points iterator
	std::vector<fitness_vector>::iterator it_p;

	// volume iterator - used for finding the contributor using std::lower_bound
	std::vector<double>::iterator it_sums;

	unsigned int i = 0;

	// Total sum of every box
	double V = 0.0;
	for(it_p = points.begin() ; it_p != points.end() ; ++it_p) {
		V = (sums[i++] = V + base::volume_between(*it_p, r_point));
	}

	unsigned long long M = 0; // Round counter
	unsigned long long M_sum = 0; // Total number of samples over every round so far

	// Container for the random point
	fitness_vector rnd_point(dim, 0.0);

	while(true) {
		// Get the random volume in-between [0, V] range, in order to choose the box with probability sums[i] / V
		double r = m_drng() * V;

		// Find the contributor using binary search
		it_sums = std::lower_bound(sums.begin(), sums.end(), r);
		i = std::distance(sums.begin(), it_sums);

		// Sample a point inside the 'box' (r_point, points[i])
		for(unsigned int d_idx = 0 ; d_idx < dim ; ++d_idx) {
			rnd_point[d_idx] = (points[i][d_idx] + m_drng() * (r_point[d_idx] - points[i][d_idx]));
		}

		unsigned int j = 0;
		do {
			if ( M_sum >= T ) {
				return (T * V) / static_cast<double>(n * M);
			}
			j = static_cast<unsigned int>(n * m_drng());
			++M_sum;
		} while (!(base::dom_cmp(rnd_point, points[j]) == base::DOM_CMP_B_DOMINATES_A));
		++M;
	}
}

/// Exclusive method
/**
 * This algorithm does not support this method.
 */
double bf_fpras::exclusive(const unsigned int p_idx, std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	(void)p_idx;
	(void)points;
	(void)r_point;
	pagmo_throw(value_error, "This method is not supported by the bf_fpras algorithm");
}

/// Least contributor method
/**
 * This algorithm does not support this method.
 */
unsigned int bf_fpras::least_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	(void)points;
	(void)r_point;
	pagmo_throw(value_error, "This method is not supported by the bf_fpras algorithm");
}

/// Greatest contributor method
/**
 * This algorithm does not support this method.
 */
unsigned int bf_fpras::greatest_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	(void)points;
	(void)r_point;
	pagmo_throw(value_error, "This method is not supported by the bf_fpras algorithm");
}

/// Contributions method
/**
 * As of yet, this algorithm does not support this method, even in its naive form, due to a poor handling of the dominated points.
 */
std::vector<double> bf_fpras::contributions(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	(void)points;
	(void)r_point;
	pagmo_throw(value_error, "This method is not supported by the bf_fpras algorithm");
}

/// Clone method.
base_ptr bf_fpras::clone() const
{
	return base_ptr(new bf_fpras(*this));
}

/// Algorithm name
std::string bf_fpras::get_name() const
{
	return "Hypervolume algorithm based on FPRAS";
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::bf_fpras)
