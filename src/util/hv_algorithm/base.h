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

#ifndef PAGMO_UTIL_HV_ALGORITHM_BASE_H
#define PAGMO_UTIL_HV_ALGORITHM_BASE_H

#include <iostream>
#include <string>
#include <typeinfo>

#include "../../config.h"
#include "../../exceptions.h"
#include "../../serialization.h"
#include "../../types.h"

namespace pagmo { namespace util {
/// Hypervolume algorithm namespace.
/**
 * This namespace contains all the algorithms implemented for the purpose of calculating the hypervolume indicator
 */
namespace hv_algorithm {

/// Base hypervolume algorithm class.
class base;

/// Alias for shared pointer to base algorithm.
typedef boost::shared_ptr<base> base_ptr;

class __PAGMO_VISIBLE base
{
	public:
		virtual ~base();
		virtual double compute(const std::vector<fitness_vector> &, const fitness_vector &) = 0;
		virtual void verify_before_compute(const std::vector<fitness_vector> &, const fitness_vector &) = 0;
	protected:
		void assert_maximal_reference_point(const std::vector<fitness_vector> &, const fitness_vector &);
};


} } }
#endif
