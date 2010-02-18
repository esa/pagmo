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

#ifndef PAGMO_ARCHIPELAGO_H
#define PAGMO_ARCHIPELAGO_H

#include <boost/scoped_ptr.hpp>
#include <boost/thread/barrier.hpp>
#include <vector>

#include "config.h"
#include "island.h"

namespace pagmo {

/// Archipelago class.
/**
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE archipelago
{
		typedef std::vector<island>::iterator iterator;
		typedef std::vector<island>::const_iterator const_iterator;
	public:
		/// Internal container of islands.
		typedef std::vector<island> container_type;
		archipelago();
		archipelago(const archipelago &);
		archipelago &operator=(const archipelago &);
		~archipelago();
		void join() const;
	private:
		void reset_barrier();
	private:
		// Container of islands.
		container_type				m_container;
		// A barrier used to synchronise the start time of all islands.
		boost::scoped_ptr<boost::barrier>	m_island_sync_point;
};

}

#endif
