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

/// hypervolume class.
/**
 * This class contains all procedures that are later accessed by population class when computing hypervolume using various methods
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE hypervolume
{
	public:
		hypervolume(boost::shared_ptr<population>);
		hypervolume(const std::vector<fitness_vector> &);
		double compute(const fitness_vector &, hv_algorithm::base_ptr);

	private:
		void verify_after_construct();
		void verify_before_compute(const fitness_vector &, hv_algorithm::base_ptr);
		boost::shared_ptr<population> m_pop;
		std::vector<fitness_vector> m_points;
		fitness_vector::size_type m_f_dim;
};

} }

#endif
