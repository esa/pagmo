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

#ifndef PAGMO_UTIL_HV_ALGORITHM_FPL_H
#define PAGMO_UTIL_HV_ALGORITHM_FPL_H

#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

#include "base.h"
#include "../hypervolume.h"

#include "fpl_cpp_original/hv.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// fpl hypervolume algorithm
/**
 * This is the class containing the interface to the the hypervolume computation package FPL (version 2.0).
 *
 * The original code was altered in following ways:
 * - Line 496 ("#define VARIANT 4") was added in order to make the file self-sustainable without external makefiles. This is also the default value according to the original FPL pacage makefile.
 *
 * @see Andreia P. Guerreiro, Carlos M. Fonseca, Michael T. Emmerich, "A Fast Dimension-Sweep Algorithm for the Hypervolume Indicator in Four Dimensions", CCCG 2012, Charlottetown, P.E.I., August 8–10, 2012.
 * @see C. M. Fonseca, L. Paquete, and M. Lopez-Ibanez. "An improved dimension-sweep algorithm for the hypervolume indicator". In IEEE Congress on Evolutionary Computation, pages 1157-1163, Vancouver, Canada, July 2006.
 * @see Nicola Beume, Carlos M. Fonseca, Manuel López-Ibáñez, Luís Paquete, and J. Vahrenhold. "On the complexity of computing the hypervolume indicator". IEEE Transactions on Evolutionary Computation, 13(5):1075-1082, 2009.
 *
 * @author (original implementation) C. M. Fonseca, L. Paquete, M. Lopez-Ibanez and A. P. Guerreiro
 * @author (C++ wrapper) Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE fpl : public base
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

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::fpl)

#endif
