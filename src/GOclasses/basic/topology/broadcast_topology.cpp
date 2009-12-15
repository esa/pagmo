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

// 06/02/2009: Initial version by Marek Ruciński.

#include <utility>
#include <vector>

#include "broadcast_topology.h"

broadcast_topology::broadcast_topology():graph_topology(),m_first(0) {}

broadcast_topology::broadcast_topology(const broadcast_topology &b)
		:graph_topology(b), m_first(b.m_first) { }

broadcast_topology &broadcast_topology::operator=(const broadcast_topology &)
{
	pagmo_assert(false);
	return *this;
}

void broadcast_topology::push_back(const size_t& id)
{
	// Store frequently-used variables.
	const size_t t_size = get_number_of_nodes();

	// Add the new node
	add_node(id);

	// If this is the first node, store the id
	if (t_size == 0) {
		m_first = id;
	} else { // otherwise, connect to the first node
		add_edge(m_first, id);
		add_edge(id, m_first);
	}
}
