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

#ifndef PAGMO_ISLAND_STORAGE_H
#define PAGMO_ISLAND_STORAGE_H

#include <boost/tuple/tuple.hpp>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "algorithm/base.h"
#include "problem/base.h"
#include "rng.h"
#include "types.h"

namespace pagmo
{

/// Island storage class.
/**
 * This class holds data for the pagmo::island class, hiding the data members as private and providing a set of protected
 * methods for use by pagmo::island's friends.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE island_storage
{
	public:
		/// Individuals stored in the island are populations of tuples of decision vector, velocity vector, current fitness vector and best fitness vector.
		typedef boost::tuple<decision_vector,decision_vector,fitness_vector,fitness_vector> individual_type;
		/// Champion type.
		/**
		 * A champion is the best individual that ever lived on the island. It is defined by a decision vector and a fitness vector.
		 */
		typedef boost::tuple<decision_vector,fitness_vector> champion_type;
		/// Alias for population type.
		typedef std::vector<individual_type> population_type;
		/// Alias for island size type.
		typedef population_type::size_type size_type;
		island_storage(const problem::base &, const algorithm::base &, int n = 0);
		island_storage &operator=(const island_storage &);
	protected:
		island_storage();
		const problem::base &prob() const;
		const algorithm::base &algo() const;
		const population_type &pop() const;
		const champion_type &champion() const;
	private:
		// Data members.
		problem::base_ptr	m_prob;
		algorithm::base_ptr	m_algo;
		// Container of individuals.
		population_type		m_pop;
		// Island champion.
		champion_type		m_champion;
		// Double precision random number generator.
		mutable	rng_double	m_drng;
};

}

#endif
