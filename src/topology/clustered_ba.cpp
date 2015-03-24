/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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
#include <boost/random/uniform_int.hpp>
#include <cstddef>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>


#include "../exceptions.h"
#include "../rng.h"
#include "clustered_ba.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Constructor from kernel size and number of edges.
/**
 * Build a BA network with clustering mechanism applied. The initial kernel has size m0 and the elements
 * being inserted after the construction of the kernel is completed are connected randomly to a maximum of m nodes. These nodes
 * are then connected to eachother with probability p.
 * Will fail if m0 < 2, if m < 1 or if m > m0 or if p < 0 or p > 1
 *
 * @param[in] m0 size of the kernel
 * @param[in] m number of random connections to be established when a new node is added.
 * @param[in] p probability that a connection is established between two nodes that are adjacent to a new node.
 */
clustered_ba::clustered_ba(int m0, int m, double p):
        m_m0(boost::numeric_cast<std::size_t>(m0)),m_m(boost::numeric_cast<std::size_t>(m)), m_p(double(p)),
	m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>())
{
	if (m0 < 2 || m < 1 || m > m0) {
		pagmo_throw(value_error,"the value of m and m0 must be at least 1 and 2, and m must not be greater than m0");
	}
        if(p < 0 || p > 1) {
                pagmo_throw(value_error,"the value of p must be between 0 and 1");
        }
}

base_ptr clustered_ba::clone() const
{
        return base_ptr(new clustered_ba(*this));
}

void clustered_ba::connect(const vertices_size_type &idx)
{
	pagmo_assert(get_number_of_vertices() > 0);
	const vertices_size_type prev_size = get_number_of_vertices() - 1;
	if (prev_size < m_m0) {
		// If we had not built the initial m0 nodes, do it.
		// We want to connect the newcomer island with high probability, and make sure that
		// at least one connection exists (otherwise the island stays isolated).
		// NOTE: is it worth to make it a user-tunable parameter?
                const double prob = 0.0;
		// Flag indicating if at least 1 connection was added.
		bool connection_added = false;
		// Main loop.
		for (std::pair<v_iterator,v_iterator> vertices = get_vertices(); vertices.first != vertices.second; ++vertices.first) {
			// Do not consider the new vertex itself.
			if (*vertices.first != idx) {
				if (m_drng() < prob) {
					connection_added = true;
					// Add the connections
					add_edge(*vertices.first,idx);
					add_edge(idx,*vertices.first);
				}
			}
		}
		// If no connections were established and this is not the first island being inserted,
		// establish at least one connection with a random island other than n.
		if ((!connection_added) && (prev_size != 0)) {
			// Get a random vertex index between 0 and n_vertices - 1. Keep on repeating the procedure if by
			// chance we end up on idx again.
			boost::uniform_int<vertices_size_type> uni_int(0,get_number_of_vertices() - 1);
			vertices_size_type rnd;
                        do {
				rnd = uni_int(m_urng);
			} while (rnd == idx);
			// Add connections to the random vertex.
			add_edge(rnd,idx);
			add_edge(idx,rnd);
		}
	} else {
                // Now we need to add j edges, choosing the nodes with a probability
		// proportional to their number of connections. We keep track of the
		// connection established in order to avoid connecting twice to the same
		// node.
                // j is a random integer in the range 1 to m.
                boost::uniform_int<edges_size_type> uni_int2(1,m_m);
		std::size_t i = 0;
                std::size_t j = uni_int2(m_urng);
		std::pair<v_iterator,v_iterator> vertices;
		std::pair<a_iterator,a_iterator> adj_vertices;
                while (i < j) {
                        // Let's find the current total number of edges.
                        const edges_size_type n_edges = get_number_of_edges();
                        pagmo_assert(n_edges > 0);
                        boost::uniform_int<edges_size_type> uni_int(0,n_edges - 1 - i);
                        // Here we choose a random number between 0 and n_edges - 1 - i.
                        const edges_size_type rn = uni_int(m_urng);
                        edges_size_type n = 0;
			// Iterate over all vertices and accumulate the number of edges for each of them. Stop when the accumulated number of edges is greater
			// than rn. This is equivalent to giving a chance of connection to vertex v directly proportional to the number of edges departing from v.
			// You can think of this process as selecting a random edge among all the existing edges and connecting to the vertex from which the
			// selected edge departs.
			vertices = get_vertices();
			for (; vertices.first != vertices.second; ++vertices.first) {
				// Do not consider it_n.
				if (*vertices.first != idx) {
					adj_vertices = get_adjacent_vertices(*vertices.first);
					n += boost::numeric_cast<edges_size_type>(std::distance(adj_vertices.first,adj_vertices.second));
					if (n > rn) {
						break;
					}
				}
			}
			pagmo_assert(vertices.first != vertices.second);
			// If the candidate was not already connected, then add it.
			if (!are_adjacent(idx,*vertices.first)) {
                                // Connect to nodes that are already adjacent to idx with probability p.
                                // This step increases clustering in the network.
                                adj_vertices = get_adjacent_vertices(idx);
                                for(;adj_vertices.first != adj_vertices.second; ++adj_vertices.first) {
                                    if(m_drng() < m_p && *adj_vertices.first != *vertices.first && !are_adjacent(*adj_vertices.first,*vertices.first)) {
                                        add_edge(*adj_vertices.first, *vertices.first);
                                        add_edge(*vertices.first, *adj_vertices.first);
                                    }
                                }
                                // Connect to idx
				add_edge(*vertices.first,idx);
				add_edge(idx,*vertices.first);
				++i;
			}
		}
	}
}

/// Topology-specific human readable info.
/**
 * Will return a formatted string containing the size of the kernel, the maximum number of connections for newly-inserted nodes,
 * and the connection probability.
 *
 * @return string containing the parameters of the Clustered BA model.
 */
std::string clustered_ba::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "\tm0 = " << m_m0 << '\n';
	oss << "\tm = " << m_m << '\n';
        oss << "\tp = " << m_p << '\n';
	return oss.str();
}

std::string clustered_ba::get_name() const
{
        return "Clustered Barabasi-Albert";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::topology::clustered_ba)
