/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#include <string>

#include "../exceptions.h"
#include "base.h"
#include "hypercube.h"

namespace pagmo { namespace topology {

/// Default constructor.
hypercube::hypercube():base() {}

base_ptr hypercube::clone() const
{
	return base_ptr(new hypercube(*this));
}

/**
 * It's a kind of magic... an incremental algorithm for constructing a hypercube by Marek.
 *
 * The algorithm takes advantage of the observation that a hypercube has a recursive structure:
 * a hypercube of dimension d (for d > 0) is composed of two hypercubes of dimension (d-1)
 * connected in paralell. In consequence, every node in an d-dimensional hypercube belongs either to the
 * "first" or to the "second" sub-hypercube of dimension (d-1).
 * For every node in an d-dimensional hypercube one can identify a series of
 * hypercubes of decreasing dimensions, d-1, d-2, ... 1 to which the node belongs. Each
 * node belongs to exactly one hypercube for each of these dimensions.
 * So when adding a new node to a partial hypercube of dimension d, for each dimension d, d-1, ... 1 one has to:
 * <ul><li>Determine to which (first or second) sub-hypercube of dimension (d-1) the node belongs</li>
 * <li>If the node belongs to the second hypercube, one has to add a connection to the coresponding node in the first hypercube</li></ul>
 */
void hypercube::connect(const vertices_size_type &n)
{
	// Store frequently-used variables.
	const vertices_size_type t_size = get_number_of_vertices();
	pagmo_assert(t_size != 0);

    // This traces the size of the current hypercube. We start from dimension 1 and increase by one.
	vertices_size_type cube_size = 2;
	// This in turn traces the position offset between the "first" and "second" sub-hypercube. This is actually always equal to cube_size / 2
	vertices_size_type jump_size = 1;

	// We consider all hypercubes of increasing sizes up to the biggest one possible to obtain with the current number of nodes
	while(jump_size <= t_size) {
		// Ordinal 0-based position of the new node in the cube with cube_size nodes (i.e. a cube of dimension lg(cube_size))
		vertices_size_type index_in_cube = n % cube_size;

		// If the node belongs to the second "sub-hypercube"
		if(index_in_cube >= jump_size) {
			// Connect the new node with the corresponding node in the "first" sub-hypercube. Note: vertices are 0-based!
			add_edge(n, t_size - jump_size - 1);
			add_edge(t_size - jump_size - 1, n);
		}

		// Proceed to the larger hypercube
		cube_size <<= 1;
		jump_size <<= 1;
	}
}

std::string hypercube::get_name() const
{
	return "Hypercube";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::topology::hypercube)
