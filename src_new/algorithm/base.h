/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

#ifndef PAGMO_ALGORITHM_BASE_H
#define PAGMO_ALGORITHM_BASE_H

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <string>
#include <typeinfo>

#include "../config.h"
#include "../population.h"
#include "../rng.h"

namespace pagmo
{
/// Algorithm namespace.
/**
 * This namespace contains all the algorithms implemented in PaGMO.
 */
namespace algorithm {

/// Base algorithm class.
class __PAGMO_VISIBLE base;

/// Alias for shared pointer to base algorithm.
typedef boost::shared_ptr<base> base_ptr;

/// Base algorithm class.
/**
 * All algorithms implemented in PaGMO must derive from this base class. This base class provides each algorithm with one pagmo::rng_double
 * and one pagmo::rng_uint32 random number generators. Each algorithm must implement the base::evolve() method.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE base
{
		friend std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);
	public:
		base();
		/// Evolve method.
		/**
		 * The purpose of this method is to take a pagmo::population as input and evolve it towards the solution of the problem.
		 */
		virtual void evolve(population &) const = 0;
		/// Clone method.
		virtual base_ptr clone() const = 0;
		virtual ~base();
		std::string human_readable() const;
	protected:
		virtual std::string human_readable_extra() const;
	protected:
		/// Random number generator for double-precision floating point values.
		mutable rng_double	m_drng;
		/// Random number generator for unsigned integer values.
		mutable rng_uint32	m_urng;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);

}
}

#endif
