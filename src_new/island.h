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

// 04/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_ISLAND_H
#define PAGMO_ISLAND_H

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/tuple/tuple.hpp>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "algorithm/base.h"
#include "problem/base.h"
//#include "migration/MigrationPolicy.h"
#include "types.h"

namespace pagmo
{

// Forward declaration of archipelago class, needed to make friend.
class archipelago;

/// Island class.
class __PAGMO_VISIBLE island
{
		// Lock type alias.
		typedef boost::lock_guard<boost::mutex> lock_type;
		// Archipelago is your friend.
		friend class archipelago;
		// Algorithm is also your friend.
		friend class algorithm::base;
		// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &, const island &);
	public:
		/// Individuals stored in the island are populations of tuples of decision vector, velocity vector, current fitness vector and best fitness vector.
		typedef boost::tuple<decision_vector,decision_vector,fitness_vector,fitness_vector> individual;
		/// Alias for population type.
		typedef std::vector<individual> population;
		/// Alias for island size type.
		typedef population::size_type size_type;
		island(const island &);
		island(const problem::base &, const algorithm::base &);
		island &operator=(const island &);
		~island();
		std::string human_readable_terse() const;
		std::string human_readable() const;
		void join() const;
	private:

	private:
		// Data members.
		problem::base_ptr	m_prob;
		algorithm::base_ptr	m_algo;
		// Archipelago in which the island is inserted.
		archipelago		*m_archi;
		// Counts the total time spent by the island on evolution (in milliseconds).
		std::size_t		m_evo_time;
		// Container of individuals.
		population		m_population;
		// Mutex used to control evolution synchronisation.
		mutable boost::mutex	m_evo_mutex;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const island &);

}

#endif
