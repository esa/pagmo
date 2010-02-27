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

// 23/01/2009: Initial version by Francesco Biscani.

#include <boost/cstdint.hpp>
#include <sstream>

#include "ba_topology.h"

ba_topology::ba_topology(int m_0, int m, boost::uint32_t _seed)
		:graph_topology(),m_m_0(m_0),m_m(m),drng(_seed),seed(_seed)
{
	if (m_0 < 2 || m < 1 || m > m_0) {
		pagmo_throw(value_error,"the value of m and m_0 must be at least 1 and 2, and m must not be greater than m_0");
	}
}

/** Assignment operator for the RNG is used here on purpose - it will create the exact copy of the RNG, so the original network and it's copy will grow identically */
ba_topology::ba_topology(const ba_topology &b):graph_topology(b),m_m_0(b.m_m_0),m_m(b.m_m),seed(b.seed)
{
	drng = b.drng;
}

void ba_topology::push_back(const size_t& id)
{
	const size_t size = get_number_of_nodes();

	// Make sure the id is not already there. [MaRu] why?
	/// \todo maybe this should be controlled by the grow_topology itself?

	// Add the new node to the graph
	add_node(id);

	if (size < m_m_0) {
		// If we have not built the initial m_0 nodes, do it.
		// We want to connect the newcomer island with high probability, and make sure that
		// at least one connection exists (otherwise the island stays isolated).
		// NOTE: is it worth to make it a user-tunable parameter?
		const double prob = 0.8;

		// Flag indicating if at least 1 connection was added in the first phase.
		bool connectionAdded = false;

		const std::vector<size_t>& nodes_list = get_nodes();

		// Main loop.
		for (std::vector<size_t>::const_iterator it = nodes_list.begin(); it != nodes_list.end(); ++it) {
			if ((*it) != id) { // Do not consider the node itself!
				if (drng() < prob) {
					connectionAdded = true;
					// Add the connections
					add_edge(*it, id);
					add_edge(id, *it);
				}
			}
		}

		// If no connections were established and this is not the first island being inserted,
		// establish at least one connection with a random island.
		if ((!connectionAdded) && (size != 0)) {
			const size_t random_id = nodes_list[(size_t)(drng() * size)];
			add_edge(random_id, id);
			add_edge(id, random_id);
		}
	} else {
		// Let's find the total number of edges.
		size_t n_edges = get_number_of_edges();

		// Now we need to add m edges, choosing the nodes with a probability
		// proportional to their number of connections. We keep track of the
		// connection established in order to avoid connecting twice to the same
		// node.

		const std::vector<size_t>& nodes_list = get_nodes();

		size_t i = 0;
		while (i < m_m) {
			const size_t rn = (size_t)(drng() * n_edges);
			size_t n = 0;

			std::vector<size_t>::const_iterator it;

			for (it = nodes_list.begin(); it != nodes_list.end(); ++it) {
				if (*it != id) { // do not cosider the node itself
					n += get_neighbours_out(*it).size();
					if (rn < n) {
						break;
					}
				}
			}

			pagmo_assert(it != nodes_list.end());

			// If id_candidate was not already connected, then add it.
			if (!are_neighbours(id, *it)) {
				add_edge(*it, id);
				add_edge(id, *it);
				++i;
			}
		}
	}
}

ba_topology& ba_topology::operator=(const ba_topology&)
{
	pagmo_assert(false);
	return *this;
}

std::string ba_topology::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << m_m_0 << "_" << m_m << "_"<< seed;
	return tmp.str();
}
