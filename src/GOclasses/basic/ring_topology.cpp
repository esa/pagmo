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

#include "ring_topology.h"

ring_topology::ring_topology():graph_topology(), growing_topology(), m_first(0), m_last(0) {}

ring_topology::ring_topology(const ring_topology &r)
	:graph_topology(r), growing_topology(), m_first(0), m_last(0) {}

ring_topology &ring_topology::operator=(const ring_topology &)
{
	pagmo_assert(false);
	return *this;
}

void ring_topology::push_back(const size_t& id)
{
	// Store frequently-used variables.
	const size_t t_size = get_number_of_nodes();
	
	switch (t_size) {
		case 0:
			// If topology is empty, update the id of the first element.
			m_first = id;
			break;
		
		case 1:
		{
			const nlt_const_iterator b = lists_out_begin();
			/// \todo pagmo_assert(id != b->first);
			// Add connections to the only existing element.
			add_edge(id, b->first);
			add_edge(b->first, id);
			break;
		}
		
		case 2:
			/// \todo pagmo_assert(m_tc.find(id) == m_tc.end());
			// Add new connections
			add_edge(m_last, id);
			add_edge(id, m_last);
			add_edge(m_first, id);
			add_edge(id, m_first);		
			break;
			
		default:
			/// \todo pagmo_assert(m_tc.find(id) == m_tc.end());
			// In general we must change the back connection of the first,
			// the forward connection of the current last, and add the new last
			// with proper connections.
			remove_edge(m_last, m_first);
			remove_edge(m_first, m_last);
			add_edge(m_last, id);
			add_edge(id, m_last);
			add_edge(m_first, id);
			add_edge(id, m_first);
	}
	// Update the id of the last island.
	m_last = id;
}
