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

#ifndef PAGMO_UTIL_HV_ALGORITHM_HV2D_H
#define PAGMO_UTIL_HV_ALGORITHM_HV2D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "base.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// hv2d hypervolume algorithm class
/**
 * This is the class containing the implementation of the hypervolume algorithm for the 2-dimensional fronts.
 * This method achieves the lower bound of n*log(n) time by sorting the initial set of points and then computing the partial areas linearly.
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE hv2d : public base
{
public:
	hv2d(const bool initial_sorting = true);
	double compute(std::vector<fitness_vector> &, const fitness_vector &) const;
	double compute(double**, unsigned int n_points, double*) const;
	std::vector<double> contributions(std::vector<fitness_vector> &, const fitness_vector &) const;

	void verify_before_compute(const std::vector<fitness_vector> &, const fitness_vector &) const;
	base_ptr clone() const;
	std::string get_name() const;

private:
	// Flag stating whether the points should be sorted in the first step of the algorithm.
	const bool m_initial_sorting;

	static bool point_pairs_cmp(const std::pair<fitness_vector, unsigned int> &, const std::pair<fitness_vector, unsigned int> &);

	static bool cmp_double_2d(double*, double*);

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<bool &>(m_initial_sorting);
	}
};

} } }

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::hv2d)

#endif
