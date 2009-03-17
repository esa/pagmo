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

// 12/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_BASE_TOPOLOGY_H
#define PAGMO_BASE_TOPOLOGY_H

#include <string>
#include <typeinfo>

#include "../../../config.h"
#include <list>

/// Base class for topologies.
/**
 * This class defines the interface that all topologies must implement.
 * \todo Rename this class in order to keep to the standard conventions.
 */
class __PAGMO_VISIBLE base_topology {
	public:
		/// Virtual Destructor.
		virtual ~base_topology() { };
		
		
		// Topology interface functions.
		
		/// Get a list of island's neighbours (outbound edges).
		virtual std::list<size_t> get_neighbours_out(const size_t&) = 0;
		
		/// Get a list of island's neighbours (inbound edges).
		virtual std::list<size_t> get_neighbours_in(const size_t&) = 0;
		
		/// Check if a pair of islands is connected.
		/**
		 * The direction of the edge must be island1 -> island2
		 */
		virtual bool are_neighbours(const size_t& island1_id, const size_t& island2_id) = 0;
		
		/// Get the number of edges in the topology
		virtual size_t get_number_of_edges() = 0;
		
		/// Get the number of nodes in the topology
		virtual size_t get_number_of_nodes() = 0;
		
		
		// Utility functions.
		
		/// Create a deep copy of the object.
		/**
		 * This method should create the exact copy of the object, including any state variables.
		 * \return Exact deep copy of the object.
		 */
		virtual base_topology *clone() const = 0;
		
		/// Get the name of the object's class.
		/** \todo Rename this method to something more sensible. */
		std::string id_name() const {return typeid(*this).name();}
};

#endif
