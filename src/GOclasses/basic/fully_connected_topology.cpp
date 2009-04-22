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

fully_connected_topology::fully_connected_topology(const double &prob):base_topology(),graph_topology(prob) {}

fully_connected_topology::fully_connected_topology(const fully_connected_topology &f):
	base_topology(f),graph_topology(f) {}

fully_connected_topology &fully_connected_topology::operator=(const fully_connected_topology &)
{
	pagmo_assert(false);
	return *this;
}

void fully_connected_topology::push_back(const island &isl)
{
	// Store frequently-used variables.
	const size_t t_size = m_tc.size(), id = isl.id();
	// Iterate over all the existing island, storing their id and adding connections to
	// the new island in the process.
	std::vector<size_t> new_connections;
	new_connections.reserve(t_size);
	const tc_iterator it_f = m_tc.end();
	for (tc_iterator it = m_tc.begin(); it != it_f; ++it) {
		pagmo_assert(id != it->first);
		new_connections.push_back(it->first);
		it->second.push_back(id);
	}
	// Insert the new island with its connections.
	m_tc.insert(std::make_pair(id,new_connections));
}

void fully_connected_topology::reset()
{
	reset_hook();
}

void fully_connected_topology::pre_evolution(island &isl)
{
	pre_hook(isl);
}

void fully_connected_topology::post_evolution(island &isl)
{
	post_hook(isl);
}
