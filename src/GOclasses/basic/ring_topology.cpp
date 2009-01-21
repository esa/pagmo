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

// 13/01/2009: Initial version by Francesco Biscani.

#include <utility>

#include "../../Functions/rng/rng.h"
#include "../../exceptions.h"
#include "archipelago.h"
#include "base_topology.h"
#include "individual.h"
#include "island.h"
#include "ring_topology.h"

ring_topology::ring_topology(const double &prob):base_topology(),m_drng(static_rng_uint32()()),m_prob(prob),m_first(0),m_last(0)
{
	if (prob < 0 || prob > 1) {
		pagmo_throw(value_error, "probability must be in the [0,1] range");
	}
}

ring_topology::ring_topology(const ring_topology &r):
	base_topology(r),m_tc(),m_ic(),m_drng(static_rng_uint32()()),m_prob(r.m_prob),m_first(0),m_last(0) {}

ring_topology &ring_topology::operator=(const ring_topology &r)
{
	if (this != &r) {
		m_tc.clear();
		m_ic.clear();
		m_drng = rng_double(static_rng_uint32()());
		m_prob = r.m_prob;
		m_first = 0;
		m_last = 0;
	}
	return *this;
}

void ring_topology::push_back(const island &isl)
{
	// Store frequently-used variables.
	const size_t t_size = m_tc.size(), id = isl.id();
	switch (t_size) {
		case 0:
			// If topology is empty, insert the id with no connections and update
			// the id of the first element.
			m_tc.insert(std::make_pair(id,std::vector<size_t>()));
			m_first = id;
			break;
		case 1:
			{
			const tc_iterator b = m_tc.begin();
			pagmo_assert(id != b->first);
			// Add a connection to the only existing element.
			b->second.push_back(id);
			// Insert new element and connect it to first.
			m_tc.insert(std::make_pair(id,std::vector<size_t>(1,b->first)));
			}
			break;
		case 2:
			pagmo_assert(m_tc.find(id) == m_tc.end());
			// The first must now have a back connection with the new last and a
			// forward connection with the current last.
			m_tc[m_first].push_back(m_tc[m_first][0]);
			m_tc[m_first][0] = id;
			// The current last must be forward connected to the new last.
			m_tc[m_last].push_back(id);
			// Insert the new last: back connection with current last and forward connection to first.
			m_tc.insert(std::make_pair(id,std::vector<size_t>(1,m_last)));
			m_tc[id].push_back(m_first);
			break;
		default:
			pagmo_assert(m_tc.find(id) == m_tc.end());
			// In general we must change the back connection of the first,
			// the forward connection of the current last, and add the new last
			// with proper connections.
			m_tc[m_first][0] = id;
			m_tc[m_last][1] = id;
			m_tc.insert(std::make_pair(id,std::vector<size_t>(1,m_last)));
			m_tc[id].push_back(m_first);
	}
	// Update the id of the last island.
	m_last = id;
}

void ring_topology::pre_evolution(island &isl)
{
	lock_type lock(m_mutex);
	// We don't want to do any migration if archipelago's size is less than 2.
	if (m_tc.size() < 2) {
		return;
	}
	if (m_drng() < m_prob) {
		std::cout << "Tentative migration\n";
		// Let's look for the island inside the topology.
		const tc_iterator tc_it = m_tc.find(isl.id());
		pagmo_assert(tc_it != m_tc.end());
		// Let's store the ids of the neighbours.
		const size_t id_prev = tc_it->second[0], id_next = tc_it->second[1];
		// Look for a best individual in the map placed there by an island neighbouring isl.
		const ic_iterator ic_its[2] = {m_ic.find(id_prev), m_ic.find(id_next)};
		// Choose randomly either the previous island or the next island.
		const ic_iterator &ic_it = ic_its[(size_t)(m_drng() * 2)];
		// If the candidate island has placed in the map an individual previously, grab it and use
		// it to replace isl's worst if it is better.
		Individual &worst = get_pop(isl).worst();
		if (ic_it != m_ic.end() && ic_it->second.getFitness() < worst.getFitness()) {
			std::cout << "Migration success\n";
			worst = ic_it->second;
			// Remove the migrated individual from the list.
			m_ic.erase(ic_it);
		}
	}
}

void ring_topology::post_evolution(island &isl)
{
	lock_type lock(m_mutex);
	// Insert in the map the best individual of the island. If the island index is already present in the island,
	// then simply update it if it is better, otherwise insert it.
	const size_t id = isl.id();
	const ic_iterator it = m_ic.find(id);
	if (it == m_ic.end()) {
		m_ic.insert(std::make_pair(id,get_pop(isl).best()));
	} else {
		const Individual &best = get_pop(isl).best();
		if (best.getFitness() < it->second.getFitness()) {
			it->second = best;
		}
	}
}
