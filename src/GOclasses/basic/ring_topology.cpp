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
#include <vector>

#include "../../exceptions.h"
#include "base_topology.h"
#include "graph_topology.h"
#include "individual.h"
#include "island.h"
#include "ring_topology.h"

ring_topology::ring_topology():base_topology(),graph_topology(),m_first(0),m_last(0) {}

ring_topology::ring_topology(const ring_topology &r):
	base_topology(r),graph_topology(r),m_first(0),m_last(0) {}

ring_topology &ring_topology::operator=(const ring_topology &)
{
	pagmo_assert(false);
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
