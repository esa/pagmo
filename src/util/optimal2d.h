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

#ifndef PAGMO_UTIL_OPTIMAL2D_H
#define PAGMO_UTIL_OPTIMAL2D_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "../config.h"
#include "../serialization.h"
#include "../pagmo.h"

namespace pagmo { namespace util {

// optimal2d class.
/**
 * This is the class containing implementation of the optimal 2D algorithm for computing hypervolume.
 * This method achieves the lower bound of n*log(n) time.
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class optimal2d
{
	public:
		static double compute(const std::vector<fitness_vector> &, const fitness_vector &);
};

bool compare_fitness(const fitness_vector &, const fitness_vector &);

}}

#endif
