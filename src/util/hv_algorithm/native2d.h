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

#ifndef PAGMO_UTIL_HV_ALGORITHM_NATIVE2D_H
#define PAGMO_UTIL_HV_ALGORITHM_NATIVE2D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "base.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// Native2D hypervolume algorithm class
/**
 * This is the class containing the implementation of the native2d algorithm for computing hypervolume.
 * This method achieves the lower bound of n*log(n) time by sorting the initial set of points and then computing the partial areas linearly.
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE native2d : public base {
	public:
		native2d(const bool initial_sorting = true);

		double compute(std::vector<fitness_vector> &, const fitness_vector &);
		unsigned int least_contributor(std::vector<fitness_vector> &, const fitness_vector &);

		void verify_before_compute(const std::vector<fitness_vector> &, const fitness_vector &);
		base_ptr clone() const;
		std::string get_name() const;

	private:
		// flag stating whether the points should be sorted in the first step of the algorithm
		const bool m_initial_sorting;

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<bool &>(m_initial_sorting);
		}
};

} } }

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::native2d);

#endif
