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

#include "hypervolume.h"
#include "hv_algorithm/base.h"
#include "hv_algorithm/hv2d.h"
#include "hv_algorithm/hv3d.h"
#include "hv_algorithm/hv4d.h"
#include "hv_algorithm/wfg.h"
#include "hv_algorithm/bf_approx.h"
#include "hv_algorithm/bf_fpras.h"
#include "hv_algorithm/hoy.h"
#include "hv_algorithm/fpl.h"

namespace pagmo { namespace util {

/// Constructor from population
/**
 * Constructs a hypervolume object, where points are elicited from the referenced population object.
 *
 * @param[in] pop reference to population object from which Pareto front is computed
 * @param[in] verify flag stating whether the points should be verified after the construction. This turns off the validation for the further computation as well, use 'set_verify' flag to alter it later.
 */
hypervolume::hypervolume(boost::shared_ptr<population> pop, const bool verify) : m_copy_points(true), m_verify(verify)
	{
	m_points.resize(pop->size());
	for (population::size_type idx = 0 ; idx < pop->size() ; ++idx) {
		m_points[idx] = fitness_vector(pop->get_individual(idx).cur_f);
	}

	if (m_verify) {
		verify_after_construct();
	}
}

/// Constructor from a vector of points
/**
 * Constructs a hypervolume object from a provided set of points.
 *
 * @param[in] points vector of points for which the hypervolume is computed
 * @param[in] verify flag stating whether the points should be verified after the construction. This turns off the validation for the further computation as well, use 'set_verify' flag to alter it later.
 */
hypervolume::hypervolume(const std::vector<fitness_vector> &points, const bool verify) : m_points(points), m_copy_points(true), m_verify(verify)
{
	if (m_verify) {
		verify_after_construct();
	}
}

/// Copy constructor.
/**
 * Will perform a deep copy of hypervolume object
 *
 * @param[in] hv hypervolume object to be copied
 */
hypervolume::hypervolume(const hypervolume &hv): m_points(hv.m_points), m_copy_points(hv.m_copy_points), m_verify(hv.m_verify) { }

/// Default constructor
/**
 * Initiates hypervolume with empty set of points.
 * Used for serialization purposes.
 */
hypervolume::hypervolume() : m_copy_points(true), m_verify(true)
{
	m_points.resize(0);
}

/// Setter for 'copy_points' flag
/**
 * Sets the hypervolume as a single use object.
 * It is used in cases where we are certain that we can alter the original set of points from the hypervolume object.
 * This is useful when we don't want to make a copy of the points first, as most algorithms alter the original set.
 *
 * This may result in unexpected behaviour when used incorrectly (e.g. requesting the computation twice out of the same object)
 *
 * @param[in] copy_points boolean value stating whether the hypervolume computation may use original set
 */
void hypervolume::set_copy_points(const bool copy_points)
{
	m_copy_points = copy_points;
}

/// Getter for 'copy_points' flag
bool hypervolume::get_copy_points()
{
	return m_copy_points;
}

/// Setter for the 'verify' flag
/**
 * Turns off the verification phase.
 * By default, the hypervolume object verifies whether certain characteristics of the point set hold, such as valid dimension sizes or a reference point that suits the minimisation.
 * In order to optimize the computation when the rules above are certain, we can turn off that phase.
 *
 * This may result in unexpected behaviour when used incorrectly (e.g. requesting the computation of empty set of points)
 *
 * @param[in] verify boolean value stating whether the hypervolume computation is to be executed without verification
 */
void hypervolume::set_verify(const bool verify)
{
	m_verify = verify;
}

/// Getter for the 'verify' flag
bool hypervolume::get_verify()
{
	return m_verify;
}

/// Verify after construct method
/**
 * Verifies whether basic requirements are met for the initial set of points.
 *
 * @throws value_error if point size is empty or when the dimensions among the points differ
 */
void hypervolume::verify_after_construct() const
{
	if ( m_points.size() == 0 ) {
		pagmo_throw(value_error, "Point set cannot be empty.");
	}
	fitness_vector::size_type f_dim = m_points[0].size();
	if (f_dim <= 1) {
		pagmo_throw(value_error, "Points of dimension > 1 required.");
	}
	for (std::vector<fitness_vector>::size_type idx = 1 ; idx < m_points.size() ; ++idx) {
		if ( m_points[idx].size() != f_dim ) {
			pagmo_throw(value_error, "All point set dimensions must be equal.");
		}
	}
}

/// Verify before compute method
/**
 * Verifies whether reference point and the hypervolume method meet certain criteria.
 *
 * @param[in] r_point fitness vector describing the reference point
 *
 * @throws value_error if reference point's and point set dimension do not agree
 */
void hypervolume::verify_before_compute(const fitness_vector &r_point, hv_algorithm::base_ptr hv_algorithm) const
{
	if ( m_points[0].size() != r_point.size() ) {
		pagmo_throw(value_error, "Point set dimensions and reference point dimension must be equal.");
	}
	hv_algorithm->verify_before_compute(m_points, r_point);
}

/// Choose the best hypervolume algorithm for given task
/**
 * Returns the best method for given hypervolume computation problem.
 * As of yet, only the dimension size is taken into account.
 */
hv_algorithm::base_ptr hypervolume::get_best_compute(const fitness_vector &r_point) const
{
	unsigned int fdim = r_point.size();
	unsigned int n = m_points.size();
	if (fdim == 2) {
		return hv_algorithm::hv2d().clone();
	} else if (fdim == 3) {
		return hv_algorithm::hv3d().clone();
	} else if (fdim == 4) {
		return hv_algorithm::hv4d().clone();
	} else if (fdim == 5 && n < 80) {
		return hv_algorithm::fpl().clone();
	} else {
		return hv_algorithm::wfg().clone();
	}
}

hv_algorithm::base_ptr hypervolume::get_best_exclusive(const unsigned int p_idx, const fitness_vector &r_point) const
{
	(void)p_idx;
	// Exclusive contribution and compute method share the same "best" set of algorithms.
	return get_best_compute(r_point);
}

hv_algorithm::base_ptr hypervolume::get_best_contributions(const fitness_vector &r_point) const
{
	unsigned int fdim = r_point.size();
	if (fdim == 2) {
		return hv_algorithm::hv2d().clone();
	} else if (fdim == 3) {
		return hv_algorithm::hv3d().clone();
	} else {
		return hv_algorithm::wfg().clone();
	}
}

/// Compute hypervolume
/**
 * Computes hypervolume for given reference point, using given algorithm object.
 *
 * @param[in] r_point fitness vector describing the reference point
 * @param[in] hv_algorithm instance of the algorithm object used for the computation
 *
 * @return value representing the hypervolume
 */
double hypervolume::compute(const fitness_vector &r_point, hv_algorithm::base_ptr hv_algorithm) const
{
	if (m_verify) {
		verify_before_compute(r_point, hv_algorithm);
	}

	// copy the initial set of points, as the algorithm may alter its contents
	if (m_copy_points) {
		std::vector<fitness_vector> points_cpy(m_points.begin(), m_points.end());
		return hv_algorithm->compute(points_cpy, r_point);
	} else {
		return hv_algorithm->compute(const_cast<std::vector<fitness_vector> &>(m_points), r_point);
	}
}

/// Compute hypervolume
/**
 * Computes hypervolume for given reference point.
 * This method chooses the hv_algorithm dynamically.
 *
 * @param[in] r_point fitness vector describing the reference point
 *
 * @return value representing the hypervolume
 */
double hypervolume::compute(const fitness_vector &r_point) const
{
	return compute(r_point, get_best_compute(r_point));
}

/// Compute exclusive contribution
/**
 * Computes exclusive hypervolume for given indivdual.
 *
 * @param[in] p_idx index of the individual for whom we compute the exclusive contribution to the hypervolume
 * @param[in] r_point fitness vector describing the reference point
 * @param[in] hv_algorithm instance of the algorithm object used for the computation
 *
 * @return value representing the hypervolume
 */
double hypervolume::exclusive(const unsigned int p_idx, const fitness_vector &r_point, hv_algorithm::base_ptr hv_algorithm) const
{
	if (m_verify) {
		verify_before_compute(r_point, hv_algorithm);
	}

	if (p_idx >= m_points.size()) {
		pagmo_throw(value_error, "Index of the individual is out of bounds.");

	}

	// copy the initial set of points, as the algorithm may alter its contents
	if (m_copy_points) {
		std::vector<fitness_vector> points_cpy(m_points.begin(), m_points.end());
		return hv_algorithm->exclusive(p_idx, points_cpy, r_point);
	} else {
		return hv_algorithm->exclusive(p_idx, const_cast<std::vector<fitness_vector> &>(m_points), r_point);
	}
}

/// Compute exclusive contribution
/**
 * Computes exclusive hypervolume for given indivdual.
 * This methods chooses the hv_algorithm dynamically.
 *
 * @param[in] p_idx index of the individual for whom we compute the exclusive contribution to the hypervolume
 * @param[in] r_point fitness vector describing the reference point
 *
 * @return value representing the hypervolume
 */
double hypervolume::exclusive(const unsigned int p_idx, const fitness_vector &r_point) const
{
	return exclusive(p_idx, r_point, get_best_exclusive(p_idx, r_point));
}

/// Find the least contributing individual
/**
 * Establishes the individual contributing the least to the total hypervolume.
 *
 * @param[in] r_point fitness vector describing the reference point
 * @param[in] hv_algorithm instance of the algorithm object used for the computation
 *
 * @return index of the least contributing point
 */
unsigned int hypervolume::least_contributor(const fitness_vector &r_point, hv_algorithm::base_ptr hv_algorithm) const
{
	if (m_verify) {
		verify_before_compute(r_point, hv_algorithm);
	}

	// Trivial case
	if (m_points.size() == 1) {
		return 0;
	}

	// copy the initial set of points, as the algorithm may alter its contents
	if (m_copy_points) {
		std::vector<fitness_vector> points_cpy(m_points.begin(), m_points.end());
		return hv_algorithm->least_contributor(points_cpy, r_point);
	} else {
		return hv_algorithm->least_contributor(const_cast<std::vector<fitness_vector> &>(m_points), r_point);
	}
}

/// Find the least contributing individual
/**
 * Establishes the individual contributing the least to the total hypervolume.
 * This method chooses the best performing hv_algorithm dynamically
 *
 * @param[in] r_point fitness vector describing the reference point
 *
 * @return index of the least contributing point
 */
unsigned int hypervolume::least_contributor(const fitness_vector &r_point) const
{
	return least_contributor(r_point, get_best_contributions(r_point));
}

/// Find the most contributing individual
/**
 * Establish the individual contributing the most to the total hypervolume.
 *
 * @param[in] r_point fitness vector describing the reference point
 * @param[in] hv_algorithm instance of the algorithm object used for the computation
 *
 * @return index of the most contributing point
 */
unsigned int hypervolume::greatest_contributor(const fitness_vector &r_point, hv_algorithm::base_ptr hv_algorithm) const
{
	if (m_verify) {
		verify_before_compute(r_point, hv_algorithm);
	}

	// copy the initial set of points, as the algorithm may alter its contents
	if (m_copy_points) {
		std::vector<fitness_vector> points_cpy(m_points.begin(), m_points.end());
		return hv_algorithm->greatest_contributor(points_cpy, r_point);
	} else {
		return hv_algorithm->greatest_contributor(const_cast<std::vector<fitness_vector> &>(m_points), r_point);
	}
}

/// Find the most contributing individual
/**
 * Establish the individual contributing the most to the total hypervolume.
 * This method chooses the best performing hv_algorithm dynamically
 *
 * @param[in] r_point fitness vector describing the reference point
 *
 * @return index of the most contributing point
 */
unsigned int hypervolume::greatest_contributor(const fitness_vector &r_point) const
{
	return greatest_contributor(r_point, get_best_contributions(r_point));
}

/// Contributions method
/**
 * This method returns the exclusive contribution to the hypervolume by every point.
 * The concrete algorithm can implement this feature optimally (as opposed to calling for the exclusive contributor in a loop).
 *
 * @param[in] r_point fitness vector describing the reference point
 * @param[in] hv_algorithm instance of the algorithm object used for the computation
 * @return vector of exclusive contributions by every point
 */
std::vector<double> hypervolume::contributions(const fitness_vector &r_point, const hv_algorithm::base_ptr hv_algorithm) const
{
	if (m_verify) {
		verify_before_compute(r_point, hv_algorithm);
	}

	// Trivial case
	if (m_points.size() == 1) {
		std::vector<double> c;
		c.push_back(hv_algorithm::base::volume_between(m_points[0], r_point));
		return c;
	}

	// copy the initial set of points, as the algorithm may alter its contents
	if (m_copy_points) {
		std::vector<fitness_vector> points_cpy(m_points.begin(), m_points.end());
		return hv_algorithm->contributions(points_cpy, r_point);
	} else {
		return hv_algorithm->contributions(const_cast<std::vector<fitness_vector> &>(m_points), r_point);
	}
}

/// Contributions method
/**
 * This method returns the exclusive contribution to the hypervolume by every point.
 * The concrete algorithm can implement this feature optimally (as opposed to calling for the exclusive contributor in a loop).
 * The hv_algorithm itself is chosen dynamically, so the best performing method is employed for given task.
 *
 * @param[in] r_point fitness vector describing the reference point
 * @return vector of exclusive contributions by every point
 */
std::vector<double> hypervolume::contributions(const fitness_vector &r_point) const
{
	return contributions(r_point, get_best_contributions(r_point));
}

/// Get expected numer of operations
/**
 * Returns the expected average amount of elementary operations for given front size (n) and dimension size (d).
 * This method is used by the approximated algorithms that fall back to exact computation.
 *
 * @param[in] n size of the front
 * @param[in] d dimension size
 *
 * @return expected number of operations for given n and d
 */
double hypervolume::get_expected_operations(const unsigned int n, const unsigned int d)
{
	if (d <= 3) {
		return d * n * log(n);  // hv3d
	} else if (d == 4) {
		return 4.0 * n * n;  // hv4d
	} else {
		return 0.0005 * d * pow(n, d * 0.5);  // exponential complexity
	}
}

/// Calculate the nadir point
/**
 * Calculates the nadir point, used as the reference point
 *
 * @param[in] epsilon value that is to be added to each objective to assure strict domination nadir point by each other point in a set
 *
 * @return value representing the hypervolume
 */
fitness_vector hypervolume::get_nadir_point(const double epsilon) const
{
	fitness_vector nadir_point(m_points[0].begin(), m_points[0].end());
	for (std::vector<fitness_vector>::size_type idx = 1 ; idx < m_points.size() ; ++ idx){
		for (fitness_vector::size_type f_idx = 0 ; f_idx < m_points[0].size() ; ++f_idx){
			// assuming minimization problem, thus maximum value by each dimension is taken
			nadir_point[f_idx] = std::max(nadir_point[f_idx], m_points[idx][f_idx]);
		}
	}
	for (fitness_vector::size_type f_idx = 0 ; f_idx < nadir_point.size() ; ++f_idx) {
		nadir_point[f_idx] += epsilon;
	}
	return nadir_point;
}


/// Get points
/**
 * Will return a vector containing the points as they were set up during construction of the hypervolume object.
 *
 * @return const reference to the vector containing the fitness_vectors representing the points in the hyperspace.
 */
const std::vector<fitness_vector> hypervolume::get_points() const
{
	return m_points;
}

/// Clone method.
/**
 * Returns a clone of the object instance.
 *
 * @return shared pointer to hypervolume class instance.
 */
hypervolume_ptr hypervolume::clone() const
{
	return hypervolume_ptr(new hypervolume(*this));
}

}}
