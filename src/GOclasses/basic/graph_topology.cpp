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

// 22/01/2009: Initial version by Francesco Biscani.

#include <iostream>

#include "../../Functions/rng/rng.h"
#include "../../exceptions.h"
#include "graph_topology.h"
#include "individual.h"
#include "island.h"

graph_topology::graph_topology():m_drng(static_rng_uint32()()) { }

graph_topology::graph_topology(const graph_topology &g):m_tc(g.m_tc),m_drng(static_rng_uint32()()) {}

graph_topology &graph_topology::operator=(const graph_topology &)
{
	pagmo_assert(false);
	return *this;
}

/*
void graph_topology::pre_hook(island &isl)
{
	lock_type lock(m_mutex);
	// We don't want to do any migration if graph size is less than 2.
	if (m_tc.size() < 2) {
		return;
	}
	if (m_drng() < m_prob) {
		// Let's look for the island inside the topology.
		const tc_iterator tc_it = m_tc.find(isl.id());
		// We want to make sure that the island is really there, otherwise the topology
		// is being used incorrectly.
		pagmo_assert(tc_it != m_tc.end());
		// If the island has no connections we have no candidates for migration, simply exit.
		const size_t n_conn = tc_it->second.size();
		if (n_conn == 0) {
			return;
		}
		// Let's choose randomly an island between the ones which are connected to isl.
		const size_t id_random = tc_it->second[(size_t)(n_conn * m_drng())];
		// Let's see if the randomly selected island has placed an individual in the best
		// individuals database.
		const ic_iterator ic_it = m_ic.find(id_random);
		// If the random island has placed in the database an individual previously, grab it and use
		// it to replace isl's worst (if it is better).
		if (ic_it != m_ic.end() && isl.t_substitute_worst(ic_it->second)) {
			// Remove the migrated individual from the list.
			m_ic.erase(ic_it);
		}
	}
}

void graph_topology::post_hook(island &isl)
{
	lock_type lock(m_mutex);
	// Insert in the database the best individual of the island. If the island index is already present in the island,
	// then simply update it if it is better, otherwise insert it.
	const size_t id = isl.id();
	const ic_iterator it = m_ic.find(id);
	const Individual best = isl.t_best();
	if (it == m_ic.end()) {
		m_ic.insert(std::make_pair(id,best));
	} else {
		if (best.getFitness() < it->second.getFitness()) {
			it->second = best;
		}
	}
}*/

std::ostream &operator<<(std::ostream &os, const graph_topology &g)
{
	graph_topology::lock_type lock(g.m_mutex);
	
	for (graph_topology::tc_type::const_iterator it = g.m_tc.begin(); it != g.m_tc.end(); ++it) {
		const size_t conn_size = it->second.size();
		os << it->first;
		if (conn_size > 0) {
			os << "->";
			for (size_t i = 0; i < conn_size; ++i) {
				os << it->second[i];
				if (i < conn_size - 1) {
					os << ',';
				}
			}
		}
		os << '\n';
	}
	return os;
}
