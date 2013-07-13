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
#include <vector>
#include "../population.h"
#include "hv_algorithm/base.h"

namespace pagmo { namespace util {

class hypervolume;

/// Alias for the shared pointer to a pagmo::util::hypervolume.
typedef boost::shared_ptr<hypervolume> hypervolume_ptr;

/// Hypervolume class.
/**
 * This class allows for setting up, and solving the hypervolume computation problems.
 * Construction of the problem is either done using the population object, or from a fixed set of points.
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE hypervolume
{
	public:
		hypervolume();
		hypervolume(const hypervolume &);
		hypervolume(boost::shared_ptr<population>);
		hypervolume(const std::vector<fitness_vector> &);
		double compute(const fitness_vector &, hv_algorithm::base_ptr);
		double exclusive(const unsigned int, const fitness_vector &, hv_algorithm::base_ptr);
		unsigned int least_contributor(const fitness_vector &, hv_algorithm::base_ptr);
		fitness_vector get_nadir_point(const double);

		hypervolume_ptr clone() const;
		const std::vector<fitness_vector> &get_points() const;

	private:
		std::vector<fitness_vector> m_points;
		void verify_after_construct();
		void verify_before_compute(const fitness_vector &, hv_algorithm::base_ptr);

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & m_points;
		}
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::util::hypervolume);

#endif
