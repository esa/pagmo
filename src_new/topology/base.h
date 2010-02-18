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
 * using a directed graph in which the nodes contain the positional index of the island inside the archipelago. The internal implementation
 * of the graph uses the Boost graph library.
 *
 * @see www.boost.org/doc/libs/release/libs/graph/doc/
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 */
class __PAGMO_VISIBLE base
{
		// Typedef for island position in archipelago.
		typedef std::vector<island>::size_type idx_type;
		// Bundled island property for boost graph.
		struct island_property
		{
			island_property() {}
			island_property(const idx_type &idx):m_idx(idx) {}
			idx_type m_idx;
		};
		// Useful shortcut typedefs for graph-related types. The graph is directed and bidirectional,
		// since we need access to both in and out edges.
		typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::bidirectionalS,island_property> graph_type;
		// Vertex iterator.
		typedef boost::graph_traits<graph_type>::vertex_iterator v_iterator;
		// In-edge iterator.
		typedef boost::graph_traits<graph_type>::in_edge_iterator ie_iterator;
		// Out-edge iterator.
		typedef boost::graph_traits<graph_type>::out_edge_iterator oe_iterator;
		// Vertex descriptor.
		typedef boost::graph_traits<graph_type>::vertex_descriptor v_descriptor;
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
		//virtual base_ptr clone() const = 0;
		virtual ~base();
		std::string human_readable_terse() const;
		std::string human_readable() const;
		/** @name Graph manipulation methods. */
		//@{
		void add_node(int);
		//@}
	protected:
		virtual std::string human_readable_extra() const;
	private:
		void check_idx(const idx_type &) const;
	private:
		graph_type m_graph;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);

}}

#endif
