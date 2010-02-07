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

#ifndef PAGMO_POPULATION_H
#define PAGMO_POPULATION_H

#include <boost/tuple/tuple.hpp>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "problem/base.h"
#include "rng.h"
#include "types.h"

namespace pagmo
{

// Forward declaration of island class, needed for friendship.
class __PAGMO_VISIBLE island;

/// Population class.
/**
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE population
{
	public:
		friend class __PAGMO_VISIBLE island;
		/// Individuals stored in the population are tuples of decision vector, velocity vector, current fitness vector and best fitness vector.
		typedef boost::tuple<decision_vector,decision_vector,fitness_vector,fitness_vector> individual_type;
		/// Champion type.
		/**
		 * A champion is the best individual that ever lived in the population. It is defined by a decision vector and a fitness vector.
		 */
		typedef boost::tuple<decision_vector,fitness_vector> champion_type;
		/// Population size type.
		typedef std::vector<individual_type>::size_type size_type;
		population(const problem::base &, int n = 0);
		population(const population &);
		population &operator=(const population &);
		const individual_type &get_individual(const size_type &) const;
		const problem::base &problem() const;
		const champion_type &champion() const;
		size_type size() const;
	private:
		population();
	private:
		typedef std::vector<individual_type> container_type;
		// Data members.
		// Problem.
		problem::base_ptr	m_prob;
		// Container of individuals.
		container_type		m_container;
		// Population champion.
		champion_type		m_champion;
		// Double precision random number generator.
		mutable	rng_double	m_drng;
};

}

#endif
