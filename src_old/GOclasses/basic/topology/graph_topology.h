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

#ifndef PAGMO_GRAPH_TOPOLOGY_H
#define PAGMO_GRAPH_TOPOLOGY_H

#include "base_topology.h"
#include <iostream>
#include <map>
#include <vector>
#include <set>

/// Simple graph implementation to be used by topologies.
/**
 * \todo Rename this class.
 */
class __PAGMO_VISIBLE graph_topology: public base_topology
{
	protected:
		typedef std::set<size_t> node_set_type;
		typedef std::map<size_t,std::vector<size_t> > neighbour_lists_type; ///< Topology container type abbreviation.
		typedef neighbour_lists_type::iterator nlt_iterator; ///< Topology container iterator type abbreviation.
		typedef neighbour_lists_type::const_iterator nlt_const_iterator; ///< Topology container const iterator type abbreviation.

	public:
		/// Default constructor.
		/** Creates a topology with no vertices nor edges */
		graph_topology() { };

		/// Copy constructor.
		/**
		 * Since topologies no longer contain references to islands, but raw IDs, they can safely be deeply copied.
		 */
		graph_topology(const graph_topology& gt)
				:base_topology(), nodes(gt.nodes), lists_out(gt.lists_out), lists_in(gt.lists_in) { }

		/// \see base_topology::clear
		virtual void clear();

		/// \see base_topology::get_neighbours_out
		virtual const std::vector<size_t> get_neighbours_out(const size_t&) const;

		/// \see base_topology::get_neighbours_in
		virtual const std::vector<size_t> get_neighbours_in(const size_t&) const;

		/// \see base_topology::are_neighbours
		virtual bool are_neighbours(const size_t& island1_id, const size_t& island2_id) const;

		/// \see base_topology::get_number_of_edges
		virtual size_t get_number_of_edges() const;

		/// \see base_topology::get_nodes
		virtual const std::vector<size_t> get_nodes() const;

		/// \see base_topology::get_number_of_nodes
		virtual size_t get_number_of_nodes() const;

		/// \see base_topology::clone
		virtual base_topology *clone() const {
			return new graph_topology(*this);
		}

		/// \see base_topology::id_object()
		virtual std::string id_object() const;

	protected:
		/// Add a node to the graph.
		/** \todo Document me! */
		void add_node(const size_t& id);

		/// Check if the node has been already added to te graph.
		/** \todo Document me! */
		bool contains_node(const size_t& id) const;

		/// Check if the node has been already added to te graph, throw an exception if failed.
		/** \todo Document me! */
		void check_node_present(const size_t& id) const;

		/// Check if the node hasn't been yet added to te graph, throw an exception if failed.
		/** \todo Document me! */
		void check_node_not_present(const size_t& id) const;

		/// Remove a node
		/** \todo Document me! */
		void remove_node(const size_t& id);

		/// Create an edge from island1 to island2 (by ids).
		void add_edge(const size_t& island1_id, const size_t& island2_id);

		/// Drop the edge from island1 to island2 (by ids).
		void remove_edge(const size_t& island1_id, const size_t& island2_id);

		/// Accessor to the outbound edges list begin iterator - allows iteration through all nodes. \todo I hate this... any other way to allow iteration through all nodes?
		nlt_const_iterator lists_out_begin() {
			return lists_out.begin();
		}

		/// Accessor to the outbound edges list end iterator - allows iteration through all nodes.
		nlt_const_iterator lists_out_end() {
			return lists_out.end();
		}

		/// Accessor to the inbound edges list begin iterator - allows iteration through all nodes.
		nlt_const_iterator lists_in_begin() {
			return lists_in.begin();
		}

		/// Accessor to the inbound edges list end iterator - allows iteration through all nodes.
		nlt_const_iterator lists_in_end() {
			return lists_in.end();
		}

	private:
		/// Dummy assignment operator.
		/** Assignment is not a valid operation for topologies - throws exception when called. */
		graph_topology &operator=(const graph_topology &);

		node_set_type						nodes; ///< Nodes present in the graph.
		neighbour_lists_type				lists_out; ///< Graph structure - a map of lists of outbound edges.
		neighbour_lists_type				lists_in; ///< Graph structure - a map of lists of inbound edges.
};

#endif
