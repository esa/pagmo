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

// 22/01/2009: Initial version by Francesco Biscani.

#include "graph_topology.h"
#include "../../exceptions.h"
#include <algorithm>

graph_topology &graph_topology::operator=(const graph_topology &)
{
	pagmo_assert(false);
	return *this;
}

std::list<size_t> graph_topology::get_neighbours_out(const size_t& _island_id)
{
	return std::list<size_t>(lists_out[_island_id].begin(), lists_out[_island_id].end());
}

std::list<size_t> graph_topology::get_neighbours_in(const size_t& _island_id)
{
	return std::list<size_t>(lists_in[_island_id].begin(), lists_in[_island_id].end());
}

bool graph_topology::are_neighbours(const size_t& island1_id, const size_t& island2_id)
{
	return find(lists_out[island1_id].begin(), lists_out[island1_id].end(), island2_id) != lists_out[island1_id].end();
}

size_t graph_topology::get_number_of_edges()
{
	size_t result = 0;
	for(nlt_iterator iter = lists_out.begin(); iter != lists_out.end(); ++iter) {
		result += iter->second.size();
	}
	return result;
}

size_t graph_topology::get_number_of_nodes()
{
	return lists_out.size();
}

void graph_topology::add_edge(const size_t& island1_id, const size_t& island2_id)
{
	lists_out[island1_id].push_back(island2_id);
	lists_in[island2_id].push_back(island1_id);	
}

void graph_topology::remove_edge(const size_t& island1_id, const size_t& island2_id)
{
	std::vector<size_t>::iterator position_out = find(lists_out[island1_id].begin(), lists_out[island1_id].end(), island2_id);
	std::vector<size_t>::iterator position_in = find(lists_in[island2_id].begin(), lists_in[island2_id].end(), island1_id);
	
	lists_out[island1_id].erase(position_out);
	lists_in[island2_id].erase(position_in);
}

void graph_topology::clear_all_edges()
{
	lists_out.clear();
	lists_in.clear();
}

std::ostream &operator<<(std::ostream &os, const graph_topology &g)
{
	for (graph_topology::neighbour_lists_type::const_iterator it = g.lists_out.begin(); it != g.lists_out.end(); ++it) {
		const size_t conn_size = it->second.size();
		os << it->first;
		if (conn_size > 0) {
			os << "->";
			for (size_t i = 0; i < conn_size; ++i) {
				os << it->second[i];
				if (i < conn_size - 1) {
					os << ',';
				}
			}
		}
		os << '\n';
	}
	return os;
}
