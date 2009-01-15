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

ring_topology::ring_topology(const double &prob):base_topology(),m_drng(static_rng_uint32()()),m_prob(prob)
{
	if (prob < 0 || prob > 1) {
		pagmo_throw(value_error, "probability must be in the [0,1] range");
	}
}

ring_topology::ring_topology(const ring_topology &r):
	base_topology(r),m_container(),m_drng(static_rng_uint32()()),m_prob(r.m_prob) {}

void ring_topology::pre_evolution(island *isl)
{
	lock_type lock(m_mutex);
	// Here we are sure the island has an associated archipelago because of the checks in island.cpp.
	const size_t arch_size = get_arch(isl)->size();
	// We don't want to do any migration if archipelago's size is less than 2.
	if (arch_size < 2) {
		return;
	}
	if (m_drng() < m_prob) {
		// std::cout << "Tentative migration\n";
		// Look for a best individual in the map placed there by an island neighbouring isl.
		const size_t isl_index = get_arch(isl)->island_index(isl);
		size_t index_next = (isl_index + 1) % arch_size, index_prev;
		// We need this kludge here because if we are on the first island we cannot simply
		// decrease the index to find the "previous" island: we need to point to the last island.
		if (isl_index == 0) {
			index_prev = arch_size - 1;
		} else {
			index_prev = isl_index - 1;
		}
		const iterator its[2] = {m_container.find(index_prev), m_container.find(index_next)};
		// Choose randomly either the previous island or the next island.
		const iterator &it = its[(size_t)(m_drng() * 2)];
		// If the chosen island has placed in the map an individual previously, grab it and use
		// it to replace isl's worse if it is better.
		Individual &worst = get_pop(isl).worst();
		if (it != m_container.end() && it->second.getFitness() < worst.getFitness()) {
			// std::cout << "Migration success\n";
			worst = it->second;
			// Remove the migrated individual from the list.
			m_container.erase(it);
		}
	}
}

void ring_topology::post_evolution(island *isl)
{
	lock_type lock(m_mutex);
	// Insert in the map the best individual of the island. If the island index is already present in the island,
	// then simply update it if it is better, otherwise insert it.
	const size_t isl_index = get_arch(isl)->island_index(isl);
	const iterator it = m_container.find(isl_index);
	if (it == m_container.end()) {
		m_container.insert(std::make_pair(isl_index,get_pop(isl).best()));
	} else {
		const Individual &best = get_pop(isl).best();
		if (best.getFitness() < it->second.getFitness()) {
			it->second = best;
		}
	}
}
