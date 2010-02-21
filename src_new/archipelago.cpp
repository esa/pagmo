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

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/thread/barrier.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "archipelago.h"
#include "exceptions.h"
#include "island.h"
#include "algorithm/base.h"
#include "problem/base.h"
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
 *
 * @param[in] t topology that will be associated to the archipelago.
 */
archipelago::archipelago(const topology::base &t):m_island_sync_point(),m_topology(t.clone()) {}

/// Constructor from problem, algorithm, archipelago size, island sizes and topology.
/**
 * Constructs n islands of m individuals each, with assigned problem p and algorithm a, and inserts them with push_back() into the archipelago,
 * whose topology is set to t.
 *
 * @param[in] p problem which will be assigned to all islands.
 * @param[in] a algorithm which will be assigned to all islands.
 * @param[in] n number of islands.
 * @param[in] m number of individuals on each island.
 * @param[in] t topology.
 */
archipelago::archipelago(const problem::base &p, const algorithm::base &a, int n, int m, const topology::base &t):
	m_island_sync_point(),m_topology(t.clone())
{
	for (size_type i = 0; i < boost::numeric_cast<size_type>(n); ++i) {
		push_back(island(p,a,m));
	}
}

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

/// Add an island to the archipelago.
/**
 * Both the island and the archipelago will be synchronised before any operation takes place. The island will then
 * be appended to the archipelago and connected to the existing islands via topology::base::push_back().
 *
 * Will fail if check_island() returns false on isl.
 *
 * @param[in] isl island to be added to the archipelago.
 */
void archipelago::push_back(const island &isl)
{
	// NOTE: the joins are already done in check_island().
	if (!check_island(isl)) {
		pagmo_throw(value_error,"cannot push_back() incompatible island");
	}
	m_container.push_back(isl);
	// Tell the island that it is living in an archipelago now.
	m_container.back().m_archi = this;
	// Insert the island in the topology.
	m_topology->push_back(boost::numeric_cast<int>(m_container.size() - 1));
	// Reset the barrier.
	reset_barrier();
}

/// Get the size of the archipelago.
/**
 * @return the number of islands contained in the archipelago.
 */
archipelago::size_type archipelago::get_size() const
{
	join();
	return m_container.size();
}

/// Return human readable representation of the archipelago.
/**
 * Will return a formatted string containing:
 * - the number of islands,
 * - the output of topology::base::human_readable_terse() for the topology,
 * - the output of island::human_readable_terse() for each island.
 *
 * @return formatted string containing the human readable representation of the archipelago.
 */
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

/// Return a copy of the topology.
/**
 * @return a topology::base_ptr to a clone of the internal topology.
 */
topology::base_ptr archipelago::get_topology() const
{
	join();
	return m_topology->clone();
}

/// Set topology.
/**
 * A valid topology must contain all and only the island indices of the current archipelago. I.e., if the size of
 * the archipelago is n, then the vertices list of the topology being set must contain all and only the integers between 0 and n-1.
 *
 * If this condition is satisfied, then the incoming topology will become the new archipelago topology. Otherwise, a value_error exception
 * will be raised.
 *
 * @param[in] t new topology for the archipelago.
 */
void archipelago::set_topology(const topology::base &t)
{
	join();
	if (m_container.size() != boost::numeric_cast<size_type>(t.get_number_of_vertices())) {
		pagmo_throw(value_error,"invalid topology, wrong number of vertices");
	}
	// Get the list of vertices and order it in ascending order.
	std::vector<int> vertices_list(t.get_vertices());
	std::sort(vertices_list.begin(),vertices_list.end());
	// t is a valid topology only if all and only island indices are represented in the vertices list.
	for (std::vector<int>::size_type i = 0; i < vertices_list.size(); ++i) {
		if (i != boost::numeric_cast<std::vector<int>::size_type>(vertices_list[i])) {
			pagmo_throw(value_error,"invalid topology, missing island index");
		}
	}
	// The topology is ok, assign it.
	m_topology = t.clone();
}

/// Check whether an island is compatible with the archipelago.
/**
 * Will return true if any of these conditions holds:
 * - archipelago is empty,
 * - the problem of the incoming island is compatible with the problems of the islands already
 *   present in the archipelago (via problem::base::is_compatible()).
 *
 * Otherwise, false will be returned.
 *
 * @param[in] isl island which will be checked.
 *
 * @return true if isl is compatible with the archipelago, false otherwise.
 */
bool archipelago::check_island(const island &isl) const
{
	join();
	isl.join();
	// We need to perform checks only if there are other islands in the archipelago.
	// Otherwise, any island will be allowed in.
	if (m_container.size() && (!isl.m_pop.problem().is_compatible(m_container[0].m_pop.problem()))) {
		return false;
	}
	return true;
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

/// Overload stream operator for pagmo::archipelago.
/**
 * Equivalent to printing archipelago::human_readable() to stream.
 *
 * @param[in] s stream to which the archipelago will be sent.
 * @param[in] a archipelago to be sent to stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const archipelago &a)
{
	s << a.human_readable();
	return s;
}

}
