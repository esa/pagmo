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


#include "bf_approx.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// Constructor
/**
 * Constructs an instance of the algorithm
 *
 * @param[in] use_exact boolean flag stating whether algorithm is allowed to use exact algorithms for the computation
 * @param[in] trivial_subcase_size size of the sub-front (points overlapping the bounding box) for which algorithm skips to the exact computation right away
 * @param[in] eps accuracy of the approximation
 * @param[in] delta confidence of the approximation
 * @param[in] gamma constant used for computation of delta for each of the points during the sampling
 * @param[in] delta_multiplier factor with which delta diminishes each round
 * @param[in] initial_delta_coeff initial coefficient multiplied by the delta at round 0
 * @param[in] alpha coefficicient stating how accurately current lowest contributor should be sampled
 */
bf_approx::bf_approx(const bool use_exact, const unsigned int trivial_subcase_size, const double eps, const double delta, const double delta_multiplier, const double alpha, const double initial_delta_coeff, const double gamma)
	: m_use_exact(use_exact), m_trivial_subcase_size(trivial_subcase_size), m_eps(eps), m_delta(delta), m_delta_multiplier(delta_multiplier), m_alpha(alpha), m_initial_delta_coeff(initial_delta_coeff), m_gamma(gamma) { }

double bf_approx::lc_end_condition(unsigned int idx, unsigned int LC, std::vector<double> &approx_volume, std::vector<double> &point_delta)
{
	return (approx_volume[LC] + point_delta[LC]) / (approx_volume[idx] - point_delta[idx]);
}
double bf_approx::gc_end_condition(unsigned int idx, unsigned int GC, std::vector<double> &approx_volume, std::vector<double> &point_delta)
{
	return (approx_volume[idx] + point_delta[idx]) / (approx_volume[GC] - point_delta[GC]);
}
bool bf_approx::lc_erase_condition(unsigned int idx, unsigned int LC, std::vector<double> &approx_volume, std::vector<double> &point_delta)
{
	return (approx_volume[idx] - point_delta[idx]) > (approx_volume[LC] + point_delta[LC]);
}
bool bf_approx::gc_erase_condition(unsigned int idx, unsigned int GC, std::vector<double> &approx_volume, std::vector<double> &point_delta)
{
	return (approx_volume[idx] + point_delta[idx]) < (approx_volume[GC] - point_delta[GC]);
}

/// Least contributor method
/**
 * This method establishes the individual that contributes the least to the hypervolume (approximately withing given epsilon and delta).
 *
 * @param[in] points vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point distinguished "reference point".
 *
 * @return index of the least contributing point
 */
unsigned int bf_approx::least_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	return approx_extreme_contributor(points, r_point, LEAST, base::cmp_least, lc_erase_condition, lc_end_condition);
}

/// Greatest contributor method
/**
 * This method establishes the individual that contributes the most to the hypervolume (approximately withing given epsilon and delta).
 *
 * @param[in] points vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point distinguished "reference point".
 *
 * @return index of the greatest contributing point
 */
unsigned int bf_approx::greatest_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	return approx_extreme_contributor(points, r_point, GREATEST, base::cmp_greatest, gc_erase_condition, gc_end_condition);
}

/// Approximated extreme contributor method
/**
 * Compute the extreme contributor using the approximated algorithm.
 * In order to make the original algorithm work for both the least and the greatest contributor, some portions
 * of the code had to be removed to external methods.
 * The following argument and functions are passed as arguments in the corresponding least/greatest contributor methods:
 *
 * - ec_type (argument):
 *   enum stating whether given execution aims to find the least or the greatest contributor.
 *   In either scenario, certain preprocessing steps are altered to determine it faster if possible.
 *
 * - cmp_func (function):
 *   Comparison function of two contributions, stating which of the contribution values fits our purpose (lesser or greater).
 *
 * - erase_condition (function):
 *   Determines whether current extreme contributor guarantees to exceed given candidate.
 *   In such case, the other point is removed from the racing.
 *
 * - end_condition (function):
 *   Determines whether given extreme contributor guarantees be accurate withing provided epsilon error.
 *   The return value of the function is the ratio, stating the estimated error.
 */
unsigned int bf_approx::approx_extreme_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point, extreme_contrib_type ec_type, bool (*cmp_func)(double, double),
		bool (*erase_condition)(unsigned int, unsigned int, std::vector<double> &, std::vector<double> &), double (*end_condition)(unsigned int, unsigned int, std::vector<double> &, std::vector<double> &)) const
{
	m_no_samples = std::vector<unsigned long long>(points.size(), 0);
	m_no_succ_samples = std::vector<unsigned long long>(points.size(), 0);
	m_no_ops = std::vector<unsigned long long>(points.size(), 1);
	m_point_set = std::vector<unsigned int>(points.size(), 0);
	m_box_volume = std::vector<double>(points.size(), 0.0);
	m_approx_volume = std::vector<double>(points.size(), 0.0);
	m_point_delta = std::vector<double>(points.size(), 0.0);
	m_boxes = std::vector<fitness_vector>(points.size());
	m_box_points = std::vector<std::vector<unsigned int> >(points.size());

	// precomputed log factor for the point delta computation
	const double log_factor = log (2. * points.size() * (1. + m_gamma) / (m_delta * m_gamma) );

	// round counter
	unsigned int round_no = 0;

	// round delta
	double r_delta = 0.0;

	bool stop_condition = false;

	// index of extreme contributor
	unsigned int EC = 0;

	// put every point into the set
	for(unsigned int i = 0 ; i < m_point_set.size() ; ++i) {
		m_point_set[i] = i;
	}

	// Initial computation
	// - compute bounding boxes and their hypervolume
	// - set round Delta as max of hypervolumes
	// - determine points overlapping with each bounding box
	for(std::vector<fitness_vector>::size_type idx = 0 ; idx < points.size() ; ++idx) {
		m_boxes[idx] = compute_bounding_box(points, r_point, idx);
		m_box_volume[idx] = base::volume_between(points[idx], m_boxes[idx]);
		r_delta = std::max(r_delta, m_box_volume[idx]);

		for(std::vector<fitness_vector>::size_type idx2 = 0 ; idx2 < points.size() ; ++idx2) {
			if (idx == idx2) {
				continue;
			}
			int op = point_in_box(points[idx2], points[idx], m_boxes[idx]);
			if (op == 1) {
				m_box_points[idx].push_back(idx2);
			} else if (ec_type == LEAST) {
				// Execute extra checks specifically for the least contributor
				switch(op) {
					case 2:
						// since contribution by idx is guaranteed to be 0.0 (as the point is dominated) we might as well return it as the least contributor right away
						return idx;
					case 3:
						// points at idx and idx2 are equal, each of them will contribute 0.0 hypervolume, return idx as the least contributor
						return idx;
					default:
						break;
				}
			}
		}
	}

	// decrease the initial maximum volume by a constant factor
	r_delta *= m_initial_delta_coeff;

	// Main loop
	do {
		r_delta *= m_delta_multiplier;
		++round_no;

		for(unsigned int _i = 0 ; _i < m_point_set.size() ; ++_i) {
			unsigned int idx = m_point_set[_i];
			sampling_round(points, r_delta , round_no, idx, log_factor);
		}

		// sample the extreme contributor
		sampling_round(points, m_alpha * r_delta , round_no, EC, log_factor);

		// find the new extreme contributor
		for(unsigned int _i = 0 ; _i < m_point_set.size() ; ++_i) {
			unsigned int idx = m_point_set[_i];
			if(cmp_func(m_approx_volume[idx], m_approx_volume[EC])) {
				EC = idx;
			}
		}

		// erase known non-extreme contributors
		std::vector<unsigned int>::iterator it = m_point_set.begin();
		while(it != m_point_set.end()) {
			unsigned int idx = *it;
			if ((idx != EC) && erase_condition(idx, EC, m_approx_volume, m_point_delta)) {
				it = m_point_set.erase(it);
			}
			else {
				++it;
			}
		}

		// check termination condition
		stop_condition = false;
		if (m_point_set.size() <= 1) {
			stop_condition = true;
		} else {
			stop_condition = true;
			for(unsigned int _i = 0 ; _i < m_point_set.size() ; ++_i) {
				unsigned int idx = m_point_set[_i];
				if (idx == EC) {
					continue;
				}
				double d = end_condition(idx, EC, m_approx_volume, m_point_delta);
				if( d <= 0 || d > 1 + m_eps ) {
					stop_condition = false;
					break;
				}
			}
		}
	} while (!stop_condition);

	return EC;
}

/// Performs a single round of sampling for given point at index 'idx'
void bf_approx::sampling_round(const std::vector<fitness_vector> &points, const double delta, const unsigned int round, const unsigned int idx, const double log_factor) const
{
	if (m_use_exact) {
		// if the sampling for given point was already resolved using exact method
		if (m_no_ops[idx] == 0) {
			return;
		}
		
		// if the exact computation is trivial OR when the sampling takes too long in terms of elementary operations
		if ( m_box_points[idx].size() <= m_trivial_subcase_size || static_cast<double>(m_no_ops[idx]) >= hypervolume::get_expected_operations(m_box_points[idx].size(), points[0].size()) ) {
			const std::vector<unsigned int> &bp = m_box_points[idx];
			if (bp.size() == 0) {
				m_approx_volume[idx] = m_box_volume[idx];
			} else {

				int f_dim = points[0].size();

				const fitness_vector &p = points[idx];

				std::vector<fitness_vector> sub_front(bp.size(), fitness_vector(f_dim, 0.0));

				for(unsigned int p_idx = 0 ; p_idx < sub_front.size() ; ++p_idx) {
					for(unsigned int d_idx = 0 ; d_idx < sub_front[0].size() ; ++d_idx) {
						sub_front[p_idx][d_idx] = std::max( p[d_idx], points[bp[p_idx]][d_idx] );
					}
				}

				const fitness_vector &refpoint = m_boxes[idx];
				hypervolume hv_obj = hypervolume(sub_front, false);
				hv_obj.set_copy_points(false);
				double hv = hv_obj.compute(refpoint);
				m_approx_volume[idx] = m_box_volume[idx] - hv;
			}

			m_point_delta[idx] = 0.0;
			m_no_ops[idx] = 0;

			return;
		}
	}
	
	double tmp = m_box_volume[idx] / delta;
	double required_no_samples = 0.5 * ( (1. + m_gamma) * log( round ) + log_factor ) * tmp * tmp;

	while(m_no_samples[idx] < required_no_samples) {
		++m_no_samples[idx];
		if (sample_successful(points, idx)) {
			++m_no_succ_samples[idx];
		}
	}

	m_approx_volume[idx] = static_cast<double>(m_no_succ_samples[idx]) / static_cast<double>(m_no_samples[idx]) * m_box_volume[idx];
	m_point_delta[idx] = compute_point_delta(round, idx, log_factor) * m_box_volume[idx];
}

/// samples the bounding box and returns true if it fell into the exclusive hypervolume
bool bf_approx::sample_successful(const std::vector<fitness_vector> &points, const unsigned int idx) const
{
	const fitness_vector &lb = points[idx];
	const fitness_vector &ub = m_boxes[idx];
	fitness_vector rnd_p(lb.size(), 0.0);
	for(unsigned int i = 0 ; i < lb.size() ; ++ i) {
		rnd_p[i] = lb[i] + m_drng()*(ub[i]-lb[i]);
	}

	for(unsigned int i = 0 ; i < m_box_points[idx].size() ; ++i) {

		// box_p is a point overlapping the bounding box volume
		const fitness_vector &box_p = points[m_box_points[idx][i]];

		// assume that box_p DOMINATES the random point and try to prove it otherwise below
		bool dominates = true;

		// increase the number of operations by the dimension size
		m_no_ops[idx] += box_p.size() + 1;
		for(unsigned int d_idx = 0 ; d_idx < box_p.size() ; ++d_idx) {
			if(rnd_p[d_idx] < box_p[d_idx]) { // box_p does not dominate rnd_p
				dominates=false;
				break;
			}
		}
		// if the box_p dominated the rnd_p return the sample as false
		if (dominates) {
			return false;
		}
	}
	return true;
}

/// Compute delta for given point
/**
 * Uses chernoff inequality as it was proposed in the article by Bringmann and Friedrich
 * The parameters of the method are taked from the Shark implementation of the algorithm.
 */
double bf_approx::compute_point_delta(const unsigned int round_no, const unsigned int idx, const double log_factor) const
{
	return sqrt(
			0.5 * ((1. + m_gamma) * log(static_cast<double>(round_no)) + log_factor)
			/ (static_cast<double>(m_no_samples[idx]))
		);
}

/// Determine whether point 'p' influences the volume of box (a, b)
/**
 * return 0 - box (p, R) has no overlapping volume with box (a, b)
 * return 1 - box (p, R) overlaps some volume with the box (a, b)
 * return 2 - point p dominates the point a (in which case, contribution by box (a, b) is guaranteed to be 0)
 * return 3 - point p is equal to point a (box (a, b) also contributes 0 hypervolume)
 */
int bf_approx::point_in_box(const fitness_vector &p, const fitness_vector &a, const fitness_vector &b) const
{
	int cmp_a_p = base::dom_cmp(a, p);

	// point a is equal to point p (duplicate)
	if (cmp_a_p == 3) {
		return 3;
	// point a is dominated by p (a is the least contributor)
	} else if (cmp_a_p == 1) {
		return 2;
	} else if(base::dom_cmp(b, p) == 1) {
		return 1;
	} else {
		return 0;
	}
}

/// Compute bounding box method
/* Find the MINIMAL (in terms of volume) bounding box that contains all the exclusive hypervolume contributed by point at index 'p_idx'
 * Thus obtained point 'z' and the original point 'points[p_idx]' form a box in which we perform the monte carlo sampling
 *
 * @param[in] points pareto front points
 * @param[in] r_point reference point
 * @param[in] p_idx index of point for which we compute the bounding box
 *
 * @return fitness_vector describing the opposite corner of the bounding box
 */
fitness_vector bf_approx::compute_bounding_box(const std::vector<fitness_vector> &points, const fitness_vector &r_point, const unsigned int p_idx) const
{
	// z is the opposite corner of the bounding box (reference point as a 'safe' first candidate - this is the MAXIMAL bounding box as of yet)
	fitness_vector z(r_point);

	// Below we perform a reduction to the minimal bounding box.
	// We check whether given point at 'idx' is DOMINATED (strong domination) by 'p_idx' in exactly one objective, and DOMINATING (weak domination) in the remaining ones.
	const fitness_vector& p = points[p_idx];
	int worse_dim_idx;
	unsigned int f_dim = r_point.size();
	for(std::vector<fitness_vector>::size_type idx = 0 ; idx < points.size(); ++idx) {

		worse_dim_idx = -1; // initiate the possible opposite point dimension by -1 (no candidate)

		for(fitness_vector::size_type f_idx = 0; f_idx < f_dim; ++f_idx) {
			if (points[idx][f_idx] >= p[f_idx]) { // if any point is worse by given dimension, it's the potential dimension in which we reduce the box
				if (worse_dim_idx != -1) { // if given point is already worse in any previous dimension, skip to next point as it's a bad candidate
					worse_dim_idx = -1; // set the result to "no candidate" and break
					break;
				}
				worse_dim_idx = f_idx;
			}
		}
		if (worse_dim_idx != -1){ // if given point was worse only in one dimension it's the potential candidate for the bouding box reductor
			z[worse_dim_idx] = std::min(z[worse_dim_idx], points[idx][worse_dim_idx]); // reduce the bounding box
		}
	}
	return z;
}

/// Verify before compute method
/**
 * Verifies whether given algorithm suits the requested data.
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the vector of points
 *
 * @throws value_error when trying to compute the hypervolume for the non-maximal reference point
 */
void bf_approx::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	base::assert_minimisation(points, r_point);
}

/// Compute hypervolume
/**
 * This method is overloaded to throw an exception in case a hypervolume indicator computation is requested.
 *
 * @param[in] points vector of points containing the 3-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume.
 */
double bf_approx::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	(void)points;
	(void)r_point;
	pagmo_throw(value_error, "This algorithm can just approximate extreme contributions but not the hypervolume itself.");
	return 0.0;
}

/// Clone method.
base_ptr bf_approx::clone() const
{
	return base_ptr(new bf_approx(*this));
}

/// Algorithm name
std::string bf_approx::get_name() const
{
	return "Bringmann-Friedrich approximation method";
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::bf_approx)
