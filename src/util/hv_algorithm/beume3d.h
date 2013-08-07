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

#ifndef PAGMO_UTIL_HV_ALGORITHM_BEUME3D_H
#define PAGMO_UTIL_HV_ALGORITHM_BEUME3D_H

#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>

#include "base.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// Beume3D hypervolume algorithm class
/**
 * This is the class containing the implementation of the Beume3D algorithm for computing hypervolume.
 * The implementation uses std::multiset (which is based on red-black tree data structure) as a container for the sweeping front.
 * Original implementation by Beume et. al uses AVL-tree. 
 * The difference is insiginificant as the important characteristics (maintaining order, self-balancing) of both structures and the asymptotical times (O(log n) updates) are the same.
 *
 * @see "On the Complexity of Computing the Hypervolume Indicator", Nicola Beume, Carlos M. Fonseca, Manuel Lopez-Ibanez,
 * Luis Paquete, Jan Vahrenhold. IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 13, NO. 5, OCTOBER 2009
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE beume3d : public base {
	public:

		beume3d(const beume3d &orig);
		beume3d(bool initial_sorting = true);
		double compute(std::vector<fitness_vector> &, const fitness_vector &);
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

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::beume3d);

#endif
