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

// 22/01/2009: Initial version by Francesco Biscani.

#include "graph_topology.h"
#include "../../../exceptions.h"
#include <algorithm>
#include <sstream>

graph_topology &graph_topology::operator=(const graph_topology &)
{
	pagmo_assert(false);
	return *this;
}

void graph_topology::clear()
{
	nodes.clear();
	lists_out.clear();
	lists_in.clear();
}

const std::vector<size_t> graph_topology::get_neighbours_out(const size_t& _island_id) const
{
	check_node_present(_island_id);
	return lists_out.find(_island_id)->second;
}

const std::vector<size_t> graph_topology::get_neighbours_in(const size_t& _island_id) const
{
	check_node_present(_island_id);
	return lists_in.find(_island_id)->second;
}

bool graph_topology::are_neighbours(const size_t& island1_id, const size_t& island2_id) const
{
	if (contains_node(island1_id) && contains_node(island2_id)) {
		const std::vector<size_t>& out_neighbours = lists_out.find(island1_id)->second;
		return std::find(out_neighbours.begin(), out_neighbours.end(), island2_id) != out_neighbours.end();
	} else {
		return false;
	}
}

size_t graph_topology::get_number_of_edges() const
{
	size_t result = 0;
	for (nlt_const_iterator iter = lists_out.begin(); iter != lists_out.end(); ++iter) {
		result += iter->second.size();
	}
	return result;
}

const std::vector<size_t> graph_topology::get_nodes() const
{
	return std::vector<size_t>(nodes.begin(), nodes.end());
}

size_t graph_topology::get_number_of_nodes() const
{
	return nodes.size();
}

void graph_topology::add_node(const size_t& id)
{
	check_node_not_present(id);
	nodes.insert(id);
	lists_in[id] = std::vector<size_t>();
	lists_out[id] = std::vector<size_t>();
}

bool graph_topology::contains_node(const size_t& id) const
{
	return nodes.find(id) != nodes.end();
}

void graph_topology::check_node_present(const size_t& id) const
{
	if (!contains_node(id)) {
		pagmo_throw(index_error, "Topology doesn't contain the node yet!");
	}
}

void graph_topology::check_node_not_present(const size_t& id) const
{
	if (contains_node(id)) {
		pagmo_throw(index_error, "Topology already contains the node!");
	}
}

void graph_topology::remove_node(const size_t& id)
{
	check_node_present(id);
	nodes.erase(id);
}

void graph_topology::add_edge(const size_t& island1_id, const size_t& island2_id)
{
	check_node_present(island1_id);
	check_node_present(island2_id);

	lists_out[island1_id].push_back(island2_id);
	lists_in[island2_id].push_back(island1_id);
}

void graph_topology::remove_edge(const size_t& island1_id, const size_t& island2_id)
{
	check_node_present(island1_id);
	check_node_present(island2_id);

	std::vector<size_t>::iterator position_out = std::find(lists_out[island1_id].begin(), lists_out[island1_id].end(), island2_id);
	std::vector<size_t>::iterator position_in = std::find(lists_in[island2_id].begin(), lists_in[island2_id].end(), island1_id);

	lists_out[island1_id].erase(position_out);
	lists_in[island2_id].erase(position_in);
}

std::string graph_topology::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << get_number_of_nodes() << "_" << get_number_of_nodes();
	return tmp.str();
}
