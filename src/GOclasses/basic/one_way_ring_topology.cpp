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

// 06/02/2009: Initial version by Francesco Biscani.

#include <utility>
#include <vector>

#include "../../exceptions.h"
#include "base_topology.h"
#include "graph_topology.h"
#include "individual.h"
#include "island.h"
#include "one_way_ring_topology.h"

one_way_ring_topology::one_way_ring_topology(const double &prob):base_topology(),graph_topology(prob),m_first(0),m_last(0) {}

one_way_ring_topology::one_way_ring_topology(const one_way_ring_topology &r):
	base_topology(r),graph_topology(r),m_first(0),m_last(0) {}

one_way_ring_topology &one_way_ring_topology::operator=(const one_way_ring_topology &)
{
	pagmo_assert(false);
	return *this;
}

void one_way_ring_topology::push_back(const island &isl)
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
			// Add a connection from the only existing element.
			b->second.push_back(id);
			// Insert new element and connect it to first.
			m_tc.insert(std::make_pair(id,std::vector<size_t>(1,b->first)));
			}
			break;
		default:
			pagmo_assert(m_tc.find(id) == m_tc.end());
			// The current last must be connected to the new last.
			m_tc[m_last][0] = id;
			// Insert the new last with a connection to the first.
			m_tc.insert(std::make_pair(id,std::vector<size_t>(1,m_first)));
	}
	// Update the id of the last island.
	m_last = id;
}

void one_way_ring_topology::reset()
{
	reset_hook();
}

void one_way_ring_topology::pre_evolution(island &isl)
{
	pre_hook(isl);
}

void one_way_ring_topology::post_evolution(island &isl)
{
	post_hook(isl);
}
