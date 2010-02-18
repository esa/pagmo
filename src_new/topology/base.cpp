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

#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>

#include "../exceptions.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Default constructor.
/**
 * Will build an empty topology.
 */
base::base():m_graph() {};

/// Copy constructor.
/**
 * Will deep-copy all members.
 *
 * @param[in] t topology::base to be copied.
 */
base::base(const base &t):m_graph(t.m_graph) {}

/// Assignment operator.
/**
 * Will perform a deep copy of all members.
 *
 * @param[in] t topology::base which will be copied into this.
 *
 * @return reference to this.
 */
base &base::operator=(const base &t)
{
	if (this != &t) {
		m_graph = t.m_graph;
	}
	return *this;
}

/// Trivial destructor.
/**
 * No side effects.
 */
base::~base() {}

/** @name Input/output. */
//@{

/// Return terse human readable representation.
/**
 * Will return a formatted string containing:
 * - the name of the topology, in mangled C++ form,
 * - the number of nodes.
 *
 * @return string containing terse human readable representation of the topology.
 */
std::string base::human_readable_terse() const
{
	std::ostringstream s;
	s << "Topology type:\t" << typeid(*this).name() << '\n';
	s << "Number of nodes:\t" << boost::num_vertices(m_graph) << '\n';
	s << "Number of edges:\t" << boost::num_edges(m_graph) << '\n';
	return s.str();
}

/// Return complete human readable representation.
/**
 * Will return a formatted string containing:
 * - the output of human_readable_terse(),
 * - the output of human_readable_extra(),
 * - the list of nodes and edges.
 *
 * @return string containing complete human readable representation of the topology.
 */
std::string base::human_readable() const
{
	std::ostringstream s;
	s << human_readable_terse();
	s << human_readable_extra();
	std::pair<v_iterator,v_iterator> vertices = boost::vertices(m_graph);
	for (; vertices.first != vertices.second; ++vertices.first) {
		std::cout << m_graph[*vertices.first].m_idx <<  '\n';
	}
	return s.str();
}

/// Return extra information for human readable representation.
/**
 * Return extra topology-dependent information that will be displayed when calling human_readable(). Default
 * implementation will return an empty string.
 *
 * @return string containing extra information for human readable representation.
 */
std::string base::human_readable_extra() const
{
	return std::string();
}

//@}

// Check if the island index is already present in the topology. In such case, throw an exception.
void base::check_idx(const idx_type &idx) const
{
	std::pair<v_iterator,v_iterator> vs = boost::vertices(m_graph);
	if (std::find_if(vs.first,vs.second,idx_finder(m_graph,idx)) != vs.second) {
		pagmo_throw(value_error,"node already present");
	}
}

/// Add a node.
/**
 * Add node containing island positional index idx. Will fail if idx is negative or
 * if idx is already present inside the topology.
 *
 * @param[in] idx island positional index to be inserted into the topology.
 */
void base::add_node(int idx)
{
	const idx_type i = boost::numeric_cast<idx_type>(idx);
	check_idx(i);
	//m_graph[boost::add_vertex(m_graph)].m_idx = i;
	boost::add_vertex(island_property(idx),m_graph);
}


/// Overload stream insertion operator for topology::base.
/**
 * Will print to stream the output of topology::base::human_readable().
 *
 * @param[in] s stream into which topology will be inserted.
 * @param[in] t topology to be inserted into stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const base &t)
{
	s << t.human_readable();
	return s;
}

}}
