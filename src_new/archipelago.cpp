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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/thread/barrier.hpp>
#include <iostream>
#include <sstream>
#include <string>

#include "archipelago.h"
#include "exceptions.h"
#include "island.h"
#include "topology/base.h"
#include "topology/unconnected.h"

namespace pagmo {

/// Default constructor.
/**
 * Will construct an empty archipelago with topology::unconnected topology.
 */
archipelago::archipelago():m_island_sync_point(),m_topology(new topology::unconnected()) {}

/// Constructor from topology.
/**
 * Will construct an empty archipelago with provided topology (which will be deep-copied internally).
 */
archipelago::archipelago(const topology::base &t):m_island_sync_point(),m_topology(t.clone()) {}

/// Copy constructor.
/**
 * Will synchronise a before deep-copying all its elements.
 *
 * @param[in] a archipelago to be copied.
 */
archipelago::archipelago(const archipelago &a)
{
	a.join();
	m_container = a.m_container;
	m_topology = a.m_topology->clone();
	reset_barrier();
}

/// Assignment operator.
/**
 * Will synchronise this and a before deep-copying all elements from a into this.
 *
 * @param[in] a archipelago used for assignment.
 *
 * @return reference to this.
 */
archipelago &archipelago::operator=(const archipelago &a)
{
	if (this != &a) {
		join();
		a.join();
		m_container = a.m_container;
		m_topology = a.m_topology->clone();
		reset_barrier();
	}
	return *this;
}

/// Destructor.
/**
 * Will call join() before returning. No other side effects.
 */
archipelago::~archipelago()
{
	join();
}

/// Wait until evolution on each island has terminated.
/**
 * Will call iteratively island::join() on all islands of the archipelago.
 */
void archipelago::join() const
{
	const const_iterator it_f = m_container.end();
	for (const_iterator it = m_container.begin(); it != it_f; ++it) {
		it->join();
	}
}

void archipelago::push_back(const island &isl)
{
	join();
	isl.join();
	check_island(isl);
	m_container.push_back(isl);
	// Tell the island that it is living in an archipelago now.
	m_container.back().m_archi = this;
	// Insert the island in the topology.
	m_topology->push_back(boost::numeric_cast<int>(m_container.size() - 1));
	// Reset the barrier.
	reset_barrier();
}

std::string archipelago::human_readable() const
{
	join();
	std::ostringstream oss;
	oss << "Archipelago\n";
	oss << "===========\n\n";
	oss << "Number of islands:\t" << m_container.size() << "\n\n";
	oss << m_topology->human_readable_terse() << '\n';
	for (size_type i = 0; i < m_container.size(); ++i) {
		oss << "Island index: #" << i << "\n\n";
		oss << m_container[i].human_readable_terse() << '\n';
	}
	return oss.str();
}

// Check an island for insertion in the archipelago. Throw if the island is not compatible.
// Do NOT use is the archipelago has not been joined!
void archipelago::check_island(const island &isl) const
{
	// We need to perform checks only if there are other islands in the archipelago.
	// Otherwise, any island will be allowed in.
	if (m_container.size()) {
		if (!isl.m_pop.problem().is_compatible(m_container[0].m_pop.problem())) {
			pagmo_throw(value_error,"cannot push_back() incompatible island");
		}
	}
}

// Reset island synchronisation barrier. This method is intended as a shortcut,
// do NOT use it if the archipelago has not been joined!
void archipelago::reset_barrier()
{
	if (m_container.size()) {
		m_island_sync_point.reset(new boost::barrier(boost::numeric_cast<unsigned int>(m_container.size())));
	} else {
		m_island_sync_point.reset(0);
	}
}

std::ostream &operator<<(std::ostream &s, const archipelago &a)
{
	s << a.human_readable();
	return s;
}

}
