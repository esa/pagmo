/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

// 23/01/2009: Initial version by Francesco Biscani.

#include "ba_topology.h"

ba_topology::ba_topology(int m_0, int m):
	base_topology(),graph_topology(),m_m_0(m_0),m_m(m)
{
	if (m_0 < 1 || m < 1 || m > m_0) {
		pagmo_throw(value_error,"the value of m and m_0 must be at least 1, and m must not be greater than m_0");
	}
}

ba_topology::ba_topology(const ba_topology &b):base_topology(b),graph_topology(b),m_m_0(b.m_m_0),m_m(b.m_m) {}

void ba_topology::push_back(const island &isl)
{
	const size_t size = m_tc.size(), id = isl.id();
	// Make sure the id is not already there.
	pagmo_assert(m_tc.find(id) == m_tc.end());
	if (size <= m_m_0) {
		// If we have not built the initial m_0 nodes, do it.
		// We want to connect the newcomer island with high probability, and make sure that
		// at least one connection exists (otherwise the island stays isolated).
		// NOTE: is it worth to make it a user-tunable parameter?
		const double prob = 0.8;
		// We want to keep a list of island IDs as we go so that we can shoot randomly
		// in case no connections were established in the main loop below.
		std::vector<size_t> id_list;
		id_list.reserve(size);
		// This will be the list of connections from the newcomer to the rest of the island.
		std::vector<size_t> new_connections;
		// Main loop.
		const tc_iterator it_f = m_tc.end();
		for (tc_iterator it = m_tc.begin(); it != it_f; ++it) {
			id_list.push_back(it->first);
			if (m_drng() < prob) {
				// Add the current item to the list of newcomer's connected islands.
				new_connections.push_back(it->first);
				// Add newcomer to the connections of the current item.
				it->second.push_back(id);
			}
		}
		// If no connections were established and this is not the first island being inserted,
		// establish at least one connection with a random island.
		if (new_connections.empty() && size != 0) {
			const size_t random_id = id_list[(size_t)(m_drng() * size)];
			new_connections.push_back(random_id);
			m_tc[random_id].push_back(id);
		}
		// Finally, add the item with associated connections.
		m_tc.insert(std::make_pair(id,new_connections));
	} else {
		// Let's find the total number of edges.
		size_t n_edges = 0;
		const tc_iterator it_f = m_tc.end();
		for (tc_iterator it = m_tc.begin(); it != it_f; ++it) {
			n_edges += it->second.size();
		}
		// Now we need to add m edges, choosing the nodes with a probability
		// proportional to their number of connections. We keep track of the
		// connection established in order to avoid connecting twice to the same
		// node.
		std::vector<size_t> connections;
		connections.reserve(m_m);
		size_t i = 0;
		while(i < m_m) {
			const size_t rn = (size_t)(m_drng() * n_edges);
			size_t n = 0;
			tc_iterator it = m_tc.begin();
			for (; it != it_f; ++it) {
				n += it->second.size();
				if (rn < n) {
					break;
				}
			}
			pagmo_assert(it != it_f);
			const size_t id_candidate = it->first;
			// If id_candidate was not already connected, then add it.
			if (std::find(connections.begin(),connections.end(),id_candidate) == connections.end()) {
				connections.push_back(id_candidate);
				it->second.push_back(id);
				++i;
			}
		}
		m_tc.insert(std::make_pair(id,connections));
	}
}

ba_topology& ba_topology::operator=(const ba_topology&)
{
	pagmo_assert(false);
	return *this;
}
