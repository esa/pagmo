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

#ifndef PAGMO_UTIL_HV_ALGORITHM_BF_APPROX_H
#define PAGMO_UTIL_HV_ALGORITHM_BF_APPROX_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>
#include "../../rng.h"

#include "base.h"

#include "../hypervolume.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// Bringmann-Friedrich approximation method
/**
 * This is the class containing the implementation of the Bringmann-Friedrich approximation method for the computation of the least contributor to the hypervolume.
 * Default values for the parameters of the algorithm were obtained from the shark implementation of the algorithm (http://image.diku.dk/shark/doxygen_pages/html/_least_contributor_approximator_8hpp_source.html)
 *
 * @see "Approximating the least hypervolume contributor: NP-hard in general, but fast in practice", Karl Bringmann, Tobias Friedrich.
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE bf_approx : public base
{
public:
	bf_approx(const bool use_exact = true, const unsigned int trivial_subcase_size = 1, const double eps = 1e-2, const double delta = 1e-6, const double delta_multiplier = 0.775, const double m_alpha = 0.2, const double initial_delta_coeff = 0.1, const double gamma = 0.25);
	double compute(std::vector<fitness_vector> &, const fitness_vector &) const;
	unsigned int least_contributor(std::vector<fitness_vector> &, const fitness_vector &) const;
	unsigned int greatest_contributor(std::vector<fitness_vector> &, const fitness_vector &) const;
	void verify_before_compute(const std::vector<fitness_vector> &, const fitness_vector &) const;
	base_ptr clone() const;
	std::string get_name() const;

private:
	inline double compute_point_delta(const unsigned int, const unsigned int, const double) const;
	inline fitness_vector compute_bounding_box(const std::vector<fitness_vector> &, const fitness_vector &, const unsigned int) const;
	inline int point_in_box(const fitness_vector &p, const fitness_vector &a, const fitness_vector &b) const;
	inline void sampling_round(const std::vector<fitness_vector>&, const double, const unsigned int, const unsigned int, const double) const;
	inline bool sample_successful(const std::vector<fitness_vector> &, const unsigned int) const;

	enum extreme_contrib_type {
		LEAST = 1,
		GREATEST = 2
	};

	unsigned int approx_extreme_contributor(std::vector<fitness_vector> &points, const fitness_vector &r_point, extreme_contrib_type, bool (*cmp_func)(double, double),
		bool (*erase_condition)(unsigned int, unsigned int, std::vector<double> &, std::vector<double> &), double (*end_condition)(unsigned int, unsigned int, std::vector<double> &, std::vector<double> &)) const;

	static double lc_end_condition(unsigned int, unsigned int, std::vector<double>&, std::vector<double>&);
	static double gc_end_condition(unsigned int, unsigned int, std::vector<double>&, std::vector<double>&);
	static bool lc_erase_condition(unsigned int, unsigned int, std::vector<double>&, std::vector<double>&);
	static bool gc_erase_condition(unsigned int, unsigned int, std::vector<double>&, std::vector<double>&);

	// flag stating whether BF approximation should use exact computation for some exclusive hypervolumes
	const bool m_use_exact;
	
	// if the number of points overlapping the bounding box is small enough we can just compute that exactly
	// following variable states the number of points for which we perform the optimization
	const unsigned int m_trivial_subcase_size;

	// accuracy of the approximation
	const double m_eps;

	// confidence of the approximation 
	const double m_delta;

	// multiplier of the round delta value
	const double m_delta_multiplier;

	// alpha coefficient used for pushing on the sampling of the current least contributor
	const double m_alpha;

	// initial coefficient of the delta at round 0
	const double m_initial_delta_coeff;

	// constant used for the computation of point delta
	const double m_gamma;

	mutable rng_double	m_drng;

	/**
	 * 'least_contributor' method variables section
	 *
	 * Section below contains member variables that are relevant only to the least_contributor method.
	 * They are not serialized as the members below are irrelevant outside of that scope.
	 */
	// number of elementary operations performed for each point
	mutable std::vector<unsigned long long> m_no_ops;

	// number of samples for given box
	mutable std::vector<unsigned long long> m_no_samples;

	// number of "successful" samples that fell into the exclusive hypervolume
	mutable std::vector<unsigned long long> m_no_succ_samples;

	// stores the indices of points that were not yet removed during the progress of the algorithm
	mutable std::vector<unsigned int> m_point_set;

	// exact hypervolumes of the bounding boxes
	mutable std::vector<double> m_box_volume;

	// approximated exlusive hypervolume of each point
	mutable std::vector<double> m_approx_volume;

	// deltas computed for each point using chernoff inequality component
	mutable std::vector<double> m_point_delta;

	// pair (boxes[idx], points[idx]) form a box in which monte carlo sampling is performed
	mutable std::vector<fitness_vector> m_boxes;

	// list of indices of points that overlap the bounding box of each point
	// during monte carlo sampling it suffices to check only these points when deciding whether the sampling was "successful"
	mutable std::vector<std::vector<unsigned int> > m_box_points;
	/**
	 * End of 'least_contributor' method variables section
	 */

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<bool &>(m_use_exact);
		ar & const_cast<unsigned int &>(m_trivial_subcase_size);
		ar & const_cast<double &>(m_eps);
		ar & const_cast<double &>(m_delta);
		ar & const_cast<double &>(m_delta_multiplier);
		ar & const_cast<double &>(m_alpha);
		ar & const_cast<double &>(m_initial_delta_coeff);
		ar & const_cast<double &>(m_gamma);
		ar & m_drng;
	}
};

} } }

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::bf_approx)

#endif
