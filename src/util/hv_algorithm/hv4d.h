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

#ifndef PAGMO_UTIL_HV_ALGORITHM_HV4D_H
#define PAGMO_UTIL_HV_ALGORITHM_HV4D_H

#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

#include "base.h"
#include "../hypervolume.h"

#include "hv4d_cpp_original/hv.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// hv4d hypervolume algorithm
/**
 * This is the class containing the interface to the the hypervolume indicator in four dimensions as implemented by Andreia P. Guerreiro et al.
 *
 * The original code was altered in following ways:
 * - name of the main method was changed from "hv4d" to "guerreiro_hv4d" in order to distinguish it from the name of this class.
 * - main method was altered to return 0 hypervolume BEFORE allocating any memory in case of an empty set of points.
 *
 * Note: As of yet, only the compute method is available.
 * This is due to the fact that the original implementation assumes certain property about the input, namely that the available set of points is a non-dominated one.
 * In order to provide accurate answers to any of the hypervolume features involving exclusive contribution by given point, we need to allow for the input to have dominated points.
 * Since this is not possible, we cannot rely on the hv4d algorithm in those cases - remaining methods, which perform "naive" implementation of said features using the "compute" method, were overloaded to raise an exception.
 *
 * @see Andreia P. Guerreiro, Carlos M. Fonseca, Michael T. Emmerich, "A Fast Dimension-Sweep Algorithm for the Hypervolume Indicator in Four Dimensions", CCCG 2012, Charlottetown, P.E.I., August 8–10, 2012.
 *
 * @author (original implementation) Andreia P. Guerreiro
 * @author (C++ wrapper) Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE hv4d : public base
{
public:
	double compute(std::vector<fitness_vector> &, const fitness_vector &) const;

	void verify_before_compute(const std::vector<fitness_vector> &, const fitness_vector &) const;
	base_ptr clone() const;
	std::string get_name() const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
	}
};

} } }

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::hv4d)

#endif
