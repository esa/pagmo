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

// 13/01/2009: Initial version by Marek Ruciński.

#include "hypercube_topology.h"

hypercube_topology::hypercube_topology():graph_topology() { }

hypercube_topology::hypercube_topology(const hypercube_topology &r)
	:graph_topology(r), nodes(r.nodes) { }

hypercube_topology &hypercube_topology::operator=(const hypercube_topology &)
{
	pagmo_assert(false);
	return *this;
}

void hypercube_topology::push_back(const size_t& id)
{
	// Store frequently-used variables.
	const size_t t_size = nodes.size();
	
	// Add node to the graph
	add_node(id);
	
	// Add the node to the hypercube of dimension 0
	nodes.push_back(id);
	
	// Size of the current hypercube
	size_t dimension = 1;
	size_t cube_size = 2;
	size_t jump_size = 1;
		
	while(jump_size <= t_size) {
		// Index of the node in the cube_size-dimenional cube
		size_t index_in_cube = (t_size) % cube_size;
		
		if(index_in_cube > jump_size - 1) {
			// Connect the new node with the corresponding node in the prvious hypercube of the appropriate dimension
			add_edge(id, nodes[t_size - jump_size]);
			add_edge(nodes[t_size - jump_size], id);
		}
		
		++dimension;
		cube_size <<= 1;
		jump_size <<= 1;
	}
}
