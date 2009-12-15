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

// 12/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_BASE_TOPOLOGY_H
#define PAGMO_BASE_TOPOLOGY_H

#include <string>
#include <typeinfo>

#include "../../../config.h"
#include <vector>
#include "../../../exceptions.h"

/// Base class for topologies.
/**
 * This class defines the interface that all topologies must implement.
 * \todo Rename this class in order to keep to the standard conventions.
 */
class __PAGMO_VISIBLE base_topology
{

		/// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &, const base_topology &);

	public:
		// Creation functions
		/// Create the topology incrementally
		/**
		 * This is an optional method. If your topology class is to support incremental construction, ovveride it.
		 * \param[in] island_id Identifier of the island being added.
		 */
		virtual void push_back(const size_t& island_id) {
			(void)island_id;
			pagmo_throw(type_error, "This topology class does not support incremental construction!");
		}

		/// Clear the topology.
		/**
		 * After calling this method, the topology should be brought back to it's initial state, so that all
		 * nodes and edges have to be re-created.
		 */
		virtual void clear() = 0;


		// Topology interface functions.

		/// Get a list of island's neighbours (outbound edges).
		virtual const std::vector<size_t> get_neighbours_out(const size_t&) const = 0;

		/// Get a list of island's neighbours (inbound edges).
		virtual const std::vector<size_t> get_neighbours_in(const size_t&) const = 0;

		/// Check if a pair of islands is connected.
		/**
		 * The direction of the edge must be island1 -> island2
		 */
		virtual bool are_neighbours(const size_t& island1_id, const size_t& island2_id) const = 0;

		/// Get the number of edges in the topology.
		virtual size_t get_number_of_edges() const = 0;

		/// Get the vector of all nodes in the topology.
		virtual const std::vector<size_t> get_nodes() const = 0;

		/// Get the number of nodes in the topology.
		virtual size_t get_number_of_nodes() const = 0;


		// Utility functions.

		/// Create a deep copy of the object.
		/**
		 * This method should create the exact copy of the object, including any state variables.
		 * \return Exact deep copy of the object.
		 */
		virtual base_topology *clone() const = 0;

		/// Get the name of the object's class.
		/** Exposed to Python. \todo Rename this method to something more sensible. */
		std::string id_name() const {
			return typeid(*this).name();
		}

		/// Get the name identyfing the object (<b>not</b> the class).
		/** Exposed to Python. The string should identify the object, so that instanciations of the same class with different parameters are distinguishable. */
		virtual std::string id_object() const = 0;
};

/// Stream output operator.
std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base_topology &);

#endif
