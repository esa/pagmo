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

/**
 * It's a kind of magic... an incremental algorithm for constructing a hypercube by Marek.
 *
 * The algorithm takes advantage of the observation that a hypercube has a recursive structure:
 * a hypercube of dimension n (for n > 0) is composed of two connected in paralell hypercubes
 * of dimension (n-1). In consequence, every node in an n-dimensional hypercube belongs either to the
 * "first" or to the "second" sub-hypercube of dimension (n-1).
 * Also, for every node in an n-dimensional hypercube one can identify a series of
 * hypercubes of decreasing dimensions, n-1, n-2, ... 1 to which the node belongs. Each
 * node belongs to exactly one hypercube for each of these dimensions.
 * So when adding a new node to a partial hypercube of dimension n, for each dimension n, n-1, ... 1 one has to:
 * <ul><li>Determine to which "sub-hypercube" of dimension (n-1) the node belongs</li>
 * <li>If the node belongs to the "second" hypercube, one has to add a connection to a coresponding node in the "first" hypercube</li>
 * </ul>
 */
void hypercube_topology::push_back(const size_t& id)
{
	// Store frequently-used variables.
	const size_t t_size = nodes.size();

	// Add node to the graph
	add_node(id);

	// Add the node to the local list of all nodes
	nodes.push_back(id);

	// Traces the size of the current hypercube. We start from dimension 1 and increase.
	size_t cube_size = 2;
	// Traces the position offset between the "first" and "second" hypercube. This is actually always equal to cube_size / 2
	size_t jump_size = 1;

	// Consider all hypercubes of increasing sizes up to the biggest one possible to obtain with the current number of nodes
	while (jump_size <= t_size) {
		// Ordinal position of the new node in the cube with cube_size nodes (i.e. of dimension lg(cube_size))
		size_t index_in_cube = t_size % cube_size;

		// If the node belongs to the second "sub-hypercube"
		if (index_in_cube >= jump_size) {
			// Connect the new node with the corresponding node in the "first" sub-hypercube
			add_edge(id, nodes[t_size - jump_size]);
			add_edge(nodes[t_size - jump_size], id);
		}

		// Proceed to the larger hypercube
		cube_size <<= 1;
		jump_size <<= 1;
	}
}
