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
#include <cstddef>
#include <iostream>
#include <string>

#include "config.h"
#include "algorithm/base.h"
#include "island_storage.h"
#include "problem/base.h"
//#include "migration/MigrationPolicy.h"
#include "types.h"

namespace pagmo
{

// Forward declaration of archipelago class, needed to make friend.
class __PAGMO_VISIBLE archipelago;

/// Island class.
/**
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE island: public island_storage
{
		// Lock type alias.
		typedef boost::lock_guard<boost::mutex> lock_type;
		// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &, const island &);
		// Algorithm is your friend.
		friend __PAGMO_VISIBLE class algorithm::base;
	public:
		island(const island &);
		island(const problem::base &, const algorithm::base &, int n = 0);
		island &operator=(const island &);
		~island();
		std::string human_readable_terse() const;
		std::string human_readable() const;
		void join() const;
	private:
		// Archipelago that, it not null, contains the island.
		archipelago		*m_archi;
		// Counts the total time spent by the island on evolution (in milliseconds).
		std::size_t		m_evo_time;
		// Mutex used to control evolution synchronisation.
		mutable boost::mutex	m_evo_mutex;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const island &);

}

#endif
