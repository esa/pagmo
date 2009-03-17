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

#ifndef PAGMO_INSTANT_TOPOLOGY_H
#define PAGMO_INSTANT_TOPOLOGY_H

#include <vector>
#include "../../exceptions.h"

/// Base class for instant-creation topologies.
/**
 * Some topologies are restricted to have only certain numbers of nodes. Typical example of such a topology is
 * a hypercube, which can contain only numbers of edges which are powers of 2. For such topologies
 * it is usually troublesome to implement "incremental growth".
 */
class __PAGMO_VISIBLE instant_topology {
	public:
		/// Virtual Destructor.
		virtual ~instant_topology() { }
		
		/// Create a topology with the given islands. 
		/**
		 * The size list of islands is verified against being compatible with the topology.
		 */
		void initialize(const std::vector<size_t>& ids)
		{
			if(isValidNumberOfNodes(ids.size())) {
				do_initialize(ids);
			} else {
				pagmo_throw(type_error, "Instant creation topology failed due to the wrong number of islands");
			}
		}
		
		/// Check if the given number of nodes is valid for the topology
		virtual bool isValidNumberOfNodes(const size_t& nodesCount) = 0;
		 
	protected:
		/// Perform the actual topology creation (to be implemented by subclasses).
		virtual void do_initialize(const std::vector<size_t>& ids) = 0;
};

#endif
