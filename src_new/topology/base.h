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

#ifndef PAGMO_TOPOLOGY_BASE_H
#define PAGMO_TOPOLOGY_BASE_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../config.h"
#include "../island.h"

namespace pagmo {

/// Topology namespace.
/**
 * This namespace contains all the topologies implemented in PaGMO.
 */
namespace topology {

/// Base topology class.
class __PAGMO_VISIBLE base;

/// Alias for shared pointer to base topology.
typedef boost::shared_ptr<base> base_ptr;

/// Base topology class.
/**
 * This class represents a topology connecting island objects in an archipelago
 * using a directed graph in which the vertices contain the positional index of the island inside the archipelago.
 * E.g., the first island inserted into the archipelago has index 0, the second one index 1 and so on.
 *
 * The internal implementation
 * of the graph uses the Boost graph library.
 *
 * @see http://www.boost.org/doc/libs/release/libs/graph/doc/
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 */
class __PAGMO_VISIBLE base
{
	public:
		// NOTE: here it would be better to get the type out from the archipelago,
		// but can we do that? Circular dependencies...
		/// Typedef for island's positional index in an archipelago.
		typedef std::vector<island>::size_type idx_type;
	private:
		// Bundled island property for boost graph.
		struct island_property
		{
			island_property() {}
			island_property(const idx_type &idx):m_idx(idx) {}
			idx_type m_idx;
		};
	protected:
		// Useful shortcut typedefs for graph-related types. The graph is directed and bidirectional,
		// since we need access to both in and out edges.
		/// Underlying graph type for the representation of the topology.
		/**
		 * The graph is a directed and bidirectional adjacency list whose vertices embed an island_property class
		 * containing the positional index of the island in the archipelago.
		 */
		typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::bidirectionalS,island_property> graph_type;
		/// Iterator over the vertices.
		typedef boost::graph_traits<graph_type>::vertex_iterator v_iterator;
		/// Iterator over adjacent vertices.
		typedef boost::graph_traits<graph_type>::adjacency_iterator a_iterator;
		/// Iterator over inversely adjacent vertices.
		typedef graph_type::inv_adjacency_iterator ia_iterator;
		/// Vertex descriptor.
		typedef boost::graph_traits<graph_type>::vertex_descriptor v_descriptor;
		/// Vertices size type.
		typedef graph_type::vertices_size_type vertices_size_type;
		/// Edges size type.
		typedef graph_type::edges_size_type edges_size_type;
	private:
		// Helper functor to find an island idx inside the graph.
		struct idx_finder
		{
			idx_finder(const graph_type &g, const idx_type &idx):m_g(g),m_idx(idx) {}
			bool operator()(const v_descriptor &vd) const
			{
				return (m_g[vd].m_idx == m_idx);
			}
			const graph_type 	&m_g;
			const idx_type		&m_idx;
		};
	public:
		base();
		base(const base &);
		base &operator=(const base &);
		/// Clone method.
		/**
		 * Provided that the derived topology implements properly the copy constructor, virtually all implementations of this method will
		 * look like this:
@verbatim
return base_ptr(new derived_topology(*this));
@endverbatim
		 *
		 * @return topology::base_ptr to a copy of this.
		 */
		virtual base_ptr clone() const = 0;
		virtual ~base();
		std::string human_readable_terse() const;
		std::string human_readable() const;
		/** @name High-level graph access and manipulation methods. */
		//@{
		vertices_size_type get_number_of_vertices() const;
		edges_size_type get_number_of_edges() const;
		bool contains_vertex(int) const;
		void push_back(int n);
		std::vector<int> get_vertices() const;
		//@}
	protected:
		/** @name Low-level graph access and manipulation methods. */
		//@{
		v_iterator get_it(int) const;
		void add_vertex(int);
		void remove_vertex(const v_iterator &);
		bool are_adjacent(const v_iterator &, const v_iterator &) const;
		std::pair<a_iterator,a_iterator> get_adjacent_vertices(const v_iterator &) const;
		void add_edge(const v_iterator &, const v_iterator &);
		void remove_edge(const v_iterator &, const v_iterator &);
		void remove_all_edges();
		std::pair<v_iterator,v_iterator> get_vertices_it() const;
		/// Establish connections between islands during a push_back() operation.
		/**
		 * This method will be called by push_back() after the vertex corresponing to island index n has been added to the graph.
		 * The purpose of this method is to connect the newly-added vertex to other vertices according to the properties of the topology.
		 *
		 * @param[in] n positional index of the newly-inserted island.
		 */
		virtual void connect(int n) = 0;
		//@}
		virtual std::string human_readable_extra() const;
	private:
		graph_type m_graph;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);

}}

#endif
