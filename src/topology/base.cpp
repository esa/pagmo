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
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

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

/// Get name of the topology.
/**
 * Default implementation will return the class' mangled C++ name.
 *
 * @return name of the topology.
 */
std::string base::get_name() const
{
	return typeid(*this).name();
}

/// Return terse human readable representation.
/**
 * Will return a formatted string containing:
 * - the output of get_name(),
 * - the number of vertices,
 * - the output of human_readable_extra().
 *
 * @return string containing terse human readable representation of the topology.
 */
std::string base::human_readable_terse() const
{
	std::ostringstream s;
	s << "Topology type:\t" << get_name() << '\n';
	s << "\tNumber of vertices:\t" << get_number_of_vertices() << '\n';
	s << "\tNumber of edges:\t" << get_number_of_edges() << '\n';
	s << human_readable_extra() << '\n';
	return s.str();
}

/// Return complete human readable representation.
/**
 * Will return a formatted string containing:
 * - the output of human_readable_terse(),
 * - the list of vertices and edges in a DOT-like syntax.
 *
 * @return string containing complete human readable representation of the topology.
 */
std::string base::human_readable() const
{
	std::ostringstream s;
	s << human_readable_terse();
	s << "Connections:\n\n";
	std::pair<v_iterator,v_iterator> vertices = get_vertices();
	std::pair<a_iterator,a_iterator> adj_vertices;
	for (; vertices.first != vertices.second; ++vertices.first) {
		adj_vertices = get_adjacent_vertices(*vertices.first);
		s << (*vertices.first);
		if (adj_vertices.first != adj_vertices.second) {
			s << " -> {";
			while (adj_vertices.first != adj_vertices.second) {
				s << (*adj_vertices.first);
				++adj_vertices.first;
				if (adj_vertices.first != adj_vertices.second) {
					s << ' ';
				}
			}
			s << '}';
		}
		s << '\n';
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

/// Add a vertex.
/**
 * Add a new, unconnected vertex to the topology.
 */
void base::add_vertex()
{
	boost::add_vertex(m_graph);
}

// Check that a vertex number does not overflow the number of vertices in the graph.
void base::check_vertex_index(const vertices_size_type &idx) const
{
	if (idx >= num_vertices(m_graph)) {
		pagmo_throw(value_error,"invalid vertex index");
	}
}

/// Return iterator range to adjacent vertices.
/**
 * Adjacent vertices are those connected from the interested vertex.
 *
 * @param[in] idx index of the interested vertex.
 *
 * @return iterator range over the adjacent vertices.
 */
std::pair<base::a_iterator,base::a_iterator> base::get_adjacent_vertices(const vertices_size_type &idx) const
{
	check_vertex_index(idx);
	return boost::adjacent_vertices(boost::vertex(idx,m_graph),m_graph);
}

/// Return vector of adjacent vertices.
/**
 * Adjacent vertices are those connected from the interested index.
 *
 * @param[in] idx interested index.
 *
 * @return vector of adjacent indices.
 */
std::vector<base::vertices_size_type> base::get_v_adjacent_vertices(const vertices_size_type &idx) const
{
	const std::pair<base::a_iterator,base::a_iterator> tmp = get_adjacent_vertices(idx);
	return std::vector<base::vertices_size_type>(tmp.first,tmp.second);
}

/// Return true if two vertices are adjacent.
/**
 * The direction of the edge must be n -> m. Will fail if either n or m are not in the topology.
 *
 * @param[in] n first vertex.
 * @param[in] m second vertex.
 *
 * @return true if the two vertices are connected, false otherwise.
 */
bool base::are_adjacent(const vertices_size_type &n, const vertices_size_type &m) const
{
	check_vertex_index(n);
	check_vertex_index(m);
	// Find n's adjacent vertices.
	const std::pair<a_iterator,a_iterator> a_vertices = boost::adjacent_vertices(boost::vertex(n,m_graph),m_graph);
	return std::find(a_vertices.first,a_vertices.second,boost::vertex(m,m_graph)) != a_vertices.second;
}

/// Return the number of adjacent vertices.
/**
 * Adjacent vertices are those connected from the interested vertex.
 *
 * @param[in] idx index of the interested vertex.
 *
 * @return number of adjacent vertices.
 */
base::edges_size_type base::get_num_adjacent_vertices(const vertices_size_type &idx) const
{
	const std::pair<base::a_iterator,base::a_iterator> v = get_adjacent_vertices(idx);
	return boost::numeric_cast<edges_size_type>(std::distance(v.first,v.second));
}

/// Return true if two vertices are inversely adjacent.
/**
 * The direction must be m -> n. Will fail if either n or m are not in the topology.
 *
 * @param[in] n first vertex.
 * @param[in] m second vertex.
 *
 * @return true if the two vertices are inversely connected, false otherwise.
 */
bool base::are_inv_adjacent(const vertices_size_type &n, const vertices_size_type &m) const
{
	check_vertex_index(n);
	check_vertex_index(m);
	// Find n's inversely adjacent vertices.
	const std::pair<ia_iterator,ia_iterator> ia_vertices = boost::inv_adjacent_vertices(boost::vertex(n,m_graph),m_graph);
	return std::find(ia_vertices.first,ia_vertices.second,boost::vertex(m,m_graph)) != ia_vertices.second;
}

/// Return iterator range to inversely adjacent vertices.
/**
 * Inversely adjacent vertices are those connected to the interested vertex.
 *
 * @param[in] idx index of the interested vertex.
 *
 * @return iterator range over the inversely adjacent vertices.
 */
std::pair<base::ia_iterator,base::ia_iterator> base::get_inv_adjacent_vertices(const vertices_size_type &idx) const
{
	check_vertex_index(idx);
	return boost::inv_adjacent_vertices(boost::vertex(idx,m_graph),m_graph);
}

/// Return vector of inversely adjacent vertices.
/**
 * Inversely adjacent vertices are those connected to the interested index.
 *
 * @param[in] idx index of the interested vertex.
 *
 * @return vector of inversely adjacent indices.
 */
std::vector<base::vertices_size_type> base::get_v_inv_adjacent_vertices(const vertices_size_type &idx) const
{
	const std::pair<base::ia_iterator,base::ia_iterator> tmp = get_inv_adjacent_vertices(idx);
	return std::vector<base::vertices_size_type>(tmp.first,tmp.second);
}

/// Return the number of inversely adjacent vertices.
/**
 * @return number of inversely adjacent vertices.
 */
base::edges_size_type base::get_num_inv_adjacent_vertices(const vertices_size_type &idx) const
{
	const std::pair<base::ia_iterator,base::ia_iterator> v = get_inv_adjacent_vertices(idx);
	return boost::numeric_cast<edges_size_type>(std::distance(v.first,v.second));
}

/// Add an edge.
/**
 * Add an edge connecting n to m. Will fail if are_adjacent() returns true.
 *
 * @param[in] n index of the first vertex.
 * @param[in] m index of the second vertex.
 */
void base::add_edge(const vertices_size_type &n, const vertices_size_type &m)
{
	if (are_adjacent(n,m)) {
		pagmo_throw(value_error,"cannot add edge, vertices are already connected");
	}
	const std::pair<e_descriptor,bool> result = boost::add_edge(boost::vertex(n,m_graph),boost::vertex(m,m_graph),m_graph);
	pagmo_assert(result.second);
	// Assign weight 1 to the edge.
	boost::property_map<graph_type,boost::edge_weight_t>::type w = boost::get(boost::edge_weight,m_graph);
	w[result.first] = 1;
}

/// Remove an edge.
/**
 * Remove the edge connecting n to m. Will fail if are_adjacent() returns false.
 *
 * @param[in] n index of the first vertex.
 * @param[in] m index of the second vertex.
 */
void base::remove_edge(const vertices_size_type &n, const vertices_size_type &m)
{
	if (!are_adjacent(n,m)) {
		pagmo_throw(value_error,"cannot remove edge, vertices are not connected");
	}
	boost::remove_edge(boost::vertex(n,m_graph),boost::vertex(m,m_graph),m_graph);
}

/// Remove all edges.
/**
 * Remove all connections between vertices.
 */
void base::remove_all_edges()
{
	for (std::pair<v_iterator,v_iterator> vertices = get_vertices(); vertices.first != vertices.second; ++vertices.first) {
		boost::clear_vertex(*vertices.first,m_graph);
	}
}

/// Return iterator range to vertices.
/**
 * Return a pair of iterators, the first one to the first vertex of the graph, the second one to the end of the graph.
 *
 * @return begin/end iterator range for the graph.
 */
std::pair<base::v_iterator,base::v_iterator> base::get_vertices() const
{
	return boost::vertices(m_graph);
}

/// Get number of vertices.
/**
 * @return total number of vertices in the graph.
 */
base::vertices_size_type base::get_number_of_vertices() const
{
	return boost::num_vertices(m_graph);
}

/// Get number of edges.
/**
 * @return total number of edges in the graph.
 */
base::edges_size_type base::get_number_of_edges() const
{
	return boost::num_edges(m_graph);
}

/// Calculate average path length.
/**
 * Calculate and return the average path length of the underlying graph representation using Johnson's all pairs shortest paths
 * algorithm. All edges are given equal weight 1. If a node is unconnected, its distance from any other node will be the highest
 * value representable by the C++ int type. The average path length is calculated as the mean value of the shortest paths between
 * all pairs of vertices.
 *
 * @return the average path length for the topology.
 *
 * @see http://en.wikipedia.org/wiki/Johnson's_algorithm
 * @see http://www.boost.org/doc/libs/release/libs/graph/doc/johnson_all_pairs_shortest.html
 */
double base::get_average_shortest_path_length() const
{
	// Output matrix.
	std::vector<std::vector<int> > D(boost::numeric_cast<std::vector<std::vector<int> >::size_type>(get_number_of_vertices()),
		std::vector<int>(boost::numeric_cast<std::vector<int>::size_type>(get_number_of_vertices())));
	boost::johnson_all_pairs_shortest_paths(m_graph,D);
	double retval = 0;
	for (std::vector<std::vector<int> >::size_type i = 0; i < D.size(); ++i) {
		for (std::vector<int>::size_type j = 0; j < D[i].size(); ++j) {
			retval += D[i][j];
		}
	}
	if (get_number_of_vertices() < 2) {
		return 0;
	} else {
		return retval / (static_cast<double>(get_number_of_vertices()) * (get_number_of_vertices() - 1));
	}
}

/// Push back vertex.
/**
 * This method will add a vertex and will then call connect() to establish the connections between the newly-added node
 * and the existing nodes in the graph.
 */
void base::push_back()
{
	add_vertex();
	connect(get_number_of_vertices() - 1);
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
