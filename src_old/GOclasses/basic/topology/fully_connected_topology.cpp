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

// 25/02/2009: Initial version by Francesco Biscani.

#include <utility>

#include "../../../exceptions.h"
#include "../individual.h"
#include "../island.h"
#include "base_topology.h"
#include "fully_connected_topology.h"
#include "graph_topology.h"

fully_connected_topology::fully_connected_topology():graph_topology() {}

fully_connected_topology::fully_connected_topology(const fully_connected_topology &f):graph_topology(f) {}

fully_connected_topology &fully_connected_topology::operator=(const fully_connected_topology &)
{
	pagmo_assert(false);
	return *this;
}

void fully_connected_topology::push_back(const size_t& id)
{
	// Iterate over all the existing island, storing their id and adding connections to
	// the new island in the process.
	std::vector<size_t> nodes_list = get_nodes();

	// Add the new node
	add_node(id);

	for (std::vector<size_t>::const_iterator it = nodes_list.begin(); it != nodes_list.end(); ++it) {
		add_edge(id, *it);
		add_edge(*it, id);
	}
}
