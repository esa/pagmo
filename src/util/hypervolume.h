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

#ifndef PAGMO_UTIL_HYPERVOLUME_H
#define PAGMO_UTIL_HYPERVOLUME_H

#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../pagmo.h"

#include "lebmeasure.h"
#include "optimal2d.h"

namespace pagmo { namespace util {

enum hv_method { hv_lebmeasure = 0, hv_optimal2d = 1};
 
typedef double (*hv_method_prototype)(const std::vector<fitness_vector> &, const fitness_vector &);

/// hypervolume class.
/**
 * This class contains all procedures that are later accessed by population class when computing hypervolume using various methods
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class hypervolume
{
	public:
		hypervolume(const population &);
		hypervolume(const std::vector<fitness_vector> &);
		double compute(const fitness_vector &reference_point, const hv_method);
		static double lebmeasure(const std::vector<fitness_vector> &, const fitness_vector &);
		static double optimal2d(const std::vector<fitness_vector> &, const fitness_vector &);
	private:

		void verify_after_construct();
		void verify_before_compute(const fitness_vector &, const hv_method);

		static hv_method_prototype hv_methods[];
		const population *pop;
		std::vector<fitness_vector> points;
		fitness_vector::size_type f_dim;
};

} }

#endif
