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

#include <boost/scoped_ptr.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <iostream>

#include "config.h"
#include "atomic_counters/atomic_counters.h"
#include "algorithm/base.h"
#include "problem/base.h"
//#include "migration/MigrationPolicy.h"

namespace pagmo
{

// Forward declaration of archipelago class, needed to make friend.
class archipelago;

/// Island class.
class __PAGMO_VISIBLE island
{
		/// Mutex type abbreviation.
		typedef boost::mutex mutex_type;
		/// Lock guard type abbreviation.
		typedef boost::lock_guard<mutex_type> lock_type;
		/// Make friend with archipelago.
		friend class archipelago;
		/// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &, const island &);
	public:
		static atomic_counter_size_t			id_counter; ///< Counter used to generate the island ids.

		//Class fields
		size_t										m_id; ///< Island id.
		population									m_pop; ///< Island's population.
		boost::scoped_ptr<const algorithm::base>		m_goa; ///< Island's algorithm.
		archipelago									*m_a; ///< Associated archipelago (may be null).
		size_t										m_evo_time; ///< Counts the total time spent by the island on evolution (in milliseconds).
		// Mutex used to control evolution synchronisation.
		mutable mutex_type	m_evo_mutex;
};

/// Stream output operator.
std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const island &);

}

#endif
