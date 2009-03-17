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

// 25/02/2009: Initial version by Francesco Biscani.

#include <utility>

#include "../../exceptions.h"
#include "base_topology.h"
#include "fully_connected_topology.h"
#include "graph_topology.h"
#include "individual.h"
#include "island.h"

fully_connected_topology::fully_connected_topology():graph_topology(),growing_topology() {}

fully_connected_topology::fully_connected_topology(const fully_connected_topology &f):graph_topology(f),growing_topology() {}

fully_connected_topology &fully_connected_topology::operator=(const fully_connected_topology &)
{
	pagmo_assert(false);
	return *this;
}

void fully_connected_topology::push_back(const island &isl)
{
	// Store frequently-used variables.
	const size_t id = isl.id();
	// Iterate over all the existing island, storing their id and adding connections to
	// the new island in the process.
	const nlt_const_iterator it_f = lists_out_end();
	for (nlt_const_iterator it = lists_out_begin(); it != it_f; ++it) {
		pagmo_assert(id != it->first);
		add_edge(id, it->first);
		add_edge(it->first, id);		
	}
}
