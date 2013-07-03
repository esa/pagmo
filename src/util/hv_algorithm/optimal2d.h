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

#ifndef PAGMO_UTIL_HV_ALGORITHM_OPTIMAL2D_H
#define PAGMO_UTIL_HV_ALGORITHM_OPTIMAL2D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "base.h"

namespace pagmo { namespace util { namespace hv_algorithm {

// optimal2d class
/**
 * This is the class containing the implementation of the optimal2D algorithm for computing hypervolume.
 * This method achieves the lower bound of n*log(n) time by sorting the initial set of points and then compute partial areas linearly.
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE optimal2d : public base {
	public:
		double compute(const std::vector<fitness_vector> & points, const fitness_vector & reference_point);
		void verify_before_compute(const std::vector<fitness_vector> & points, const fitness_vector & reference_point);

};

inline bool compare_fitness(const fitness_vector &, const fitness_vector &);

} } }

#endif
