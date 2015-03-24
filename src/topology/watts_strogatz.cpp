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

#include <algorithm>
#include <boost/integer_traits.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cstddef>
#include <exception>
#include <string>
#include <utility>

#include "../exceptions.h"
#include "../rng.h"
#include "base.h"
#include "watts_strogatz.h"

namespace pagmo { namespace topology {

/// Constructor from K, beta and size.
/**
 * Build a Watts-Strogatz network model with K = k, rewiring probability beta and size n. Will fail if k is less than 2 or odd
 * or if the beta parameter is outside the [0,1] range.
 *
 * @param[in] k K parameter of the model.
 * @param[in] beta probability of rewiring.
 * @param[in] n size (i.e., number of vertices) of the topology.
 */
watts_strogatz::watts_strogatz(int k, const double &beta, int n):
	m_k(boost::numeric_cast<std::size_t>(k)),m_beta(beta),m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>())
{
	if (m_k < 2 || m_k % 2) {
		pagmo_throw(value_error,"the K parameter must be even and at least 2");
	}
	const vertices_size_type size = boost::numeric_cast<vertices_size_type>(n);
	// Solve potential overflow.
	if (m_k == boost::integer_traits<std::size_t>::const_max) {
		pagmo_throw(std::overflow_error,"overflow error in Watts-Strogatz model");
	}
	if (m_beta < 0 || m_beta > 1) {
		pagmo_throw(value_error,"the beta parameter must be in the [0,1] range");
	}
	if (size) {
		// Add all vertices.
		for (vertices_size_type i = 0; i < size; ++i) {
			add_vertex();
		}
		// Do the wiring.
		rewire();
	}
}

base_ptr watts_strogatz::clone() const
{
	return base_ptr(new watts_strogatz(*this));
}

// Connect an unconnected topology into a Watts-Strogatz model.
void watts_strogatz::rewire()
{
	const vertices_size_type size = get_number_of_vertices();
	pagmo_assert(size > 0);
	// First, build the ring lattice.
	for (std::pair<v_iterator,v_iterator> vertices = get_vertices(); vertices.first != vertices.second; ++vertices.first) {
		// Forward connections.
		v_iterator tmp = vertices.first;
		for (std::size_t j = 1; j <= m_k / 2; ++j) {
			++tmp;
			// Wrap around if we went past the last element.
			if (tmp == vertices.second) {
				tmp = get_vertices().first;
			}
			add_edge(*vertices.first,*tmp);
		}
		// Backward connections.
		tmp = vertices.first;
		for (std::size_t j = 1; j <= m_k / 2; ++j) {
			// Go to the last-past-one element if we are at the first.
			if (tmp == get_vertices().first) {
				tmp = vertices.second;
			}
			--tmp;
			add_edge(*vertices.first,*tmp);
		}
	}
	// Second, do the random rewires.
	boost::uniform_int<vertices_size_type> uint(0,size - 1);
	for (std::pair<v_iterator,v_iterator> vertices = get_vertices(); vertices.first != vertices.second; ++vertices.first) {
		v_iterator tmp = vertices.first;
		++tmp;
		const edges_size_type n_adj_vertices = get_num_adjacent_vertices(*vertices.first);
		// NOTE: here things are hairy, it is possible - especially near the kernel size - that
		// no rewire can take place because the node is fully connected. We must detect this, in order
		// to avoid an endless loop. That's why in the for loop we check also the number of vertices adjacent to
		// vertices.first.
		// Iterate over the next m_k / 2 vertices, without going past the end.
		for (std::size_t j = 1; tmp != vertices.second && j <= m_k / 2 && n_adj_vertices < size - 1; ++j, ++tmp) {
			if (m_drng() < m_beta) {
				// Select a random vertex that is not vertices.first itself and that would
				// not end in a duplicate edge if connected from vertices.first.
				vertices_size_type rng;
				do {
					rng = uint(m_urng);
				} while (rng == *vertices.first || are_adjacent(*vertices.first,rng));
				pagmo_assert(!are_adjacent(rng,*vertices.first));
				// Destroy edge and rewire.
				remove_edge(*vertices.first,*tmp);
				remove_edge(*tmp,*vertices.first);
				add_edge(*vertices.first,rng);
				add_edge(rng,*vertices.first);
			}
		}
	}
}

void watts_strogatz::connect(const vertices_size_type &n)
{
	const vertices_size_type size = get_number_of_vertices();
	if (size <= m_k + 1) {
		// If we are below the initial kernel size, just do a fully connected topology.
		for (std::pair<v_iterator,v_iterator> vertices = get_vertices(); vertices.first != vertices.second; ++vertices.first) {
			// Do not connect with self.
			if (n != *vertices.first) {
				add_edge(n,*vertices.first);
				add_edge(*vertices.first,n);
			}
		}
	} else {
		pagmo_assert(size > 0);
		// We need to do a complete rewire of the model.
		remove_all_edges();
		rewire();
	}
}

std::string watts_strogatz::get_name() const
{
	return "Watts-Strogatz";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::topology::watts_strogatz)
