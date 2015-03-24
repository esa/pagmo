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

#ifndef PAGMO_UTIL_HV_ALGORITHM_BF_FPRAS_H
#define PAGMO_UTIL_HV_ALGORITHM_BF_FPRAS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include "../../rng.h"

#include "base.h"

#include "../hypervolume.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// Bringmann-Friedrich approximation method
/**
 * This class contains the implementation of the Bringmann-Friedrich approximation scheme (FPRAS), reduced to a special case of approximating the hypervolume indicator.
 *
 * @see "Approximating the volume of unions and intersections of high-dimensional geometric objects", Karl Bringmann, Tobias Friedrich.
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE bf_fpras : public base
{
public:
	bf_fpras(const double eps = 1e-2, const double delta = 1e-2);

	double compute(std::vector<fitness_vector> &, const fitness_vector &) const;

	double exclusive(const unsigned int, std::vector<fitness_vector> &, const fitness_vector &) const;
	unsigned int least_contributor(std::vector<fitness_vector> &, const fitness_vector &) const;
	unsigned int greatest_contributor(std::vector<fitness_vector> &, const fitness_vector &) const;
	std::vector<double> contributions(std::vector<fitness_vector> &, const fitness_vector &) const;

	void verify_before_compute(const std::vector<fitness_vector> &, const fitness_vector &) const;
	base_ptr clone() const;
	std::string get_name() const;

private:
	// error of the approximation
	const double m_eps;
	// probabiltiy of error
	const double m_delta;

	mutable rng_double m_drng;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<double &>(m_eps);
		ar & const_cast<double &>(m_delta);
		ar & m_drng;
	}
};

} } }

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::bf_fpras)

#endif
