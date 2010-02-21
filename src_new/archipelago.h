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
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "island.h"
#include "algorithm/base.h"
#include "problem/base.h"
#include "topology/base.h"
#include "topology/unconnected.h"

namespace pagmo {

/// Archipelago class.
/**
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE archipelago
{
	public:
		/// Internal container of islands.
		typedef std::vector<island> container_type;
		/// Archipelago size type.
		typedef container_type::size_type size_type;
	private:
		// Iterators.
		typedef container_type::iterator iterator;
		typedef container_type::const_iterator const_iterator;
	public:
		archipelago();
		archipelago(const topology::base &);
		archipelago(const problem::base &p, const algorithm::base &a, int n, int m, const topology::base &t = topology::unconnected());
		archipelago(const archipelago &);
		archipelago &operator=(const archipelago &);
		~archipelago();
		void join() const;
		void push_back(const island &);
		size_type get_size() const;
		std::string human_readable() const;
		bool check_island(const island &) const;
		topology::base_ptr get_topology() const;
		void set_topology(const topology::base &);
	private:
		void reset_barrier();
	private:
		// Container of islands.
		container_type				m_container;
		// A barrier used to synchronise the start time of all islands.
		boost::scoped_ptr<boost::barrier>	m_island_sync_point;
		// Topology.
		topology::base_ptr			m_topology;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const archipelago &);

}

#endif
