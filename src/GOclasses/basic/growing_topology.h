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

#ifndef PAGMO_GROWING_TOPOLOGY_H
#define PAGMO_GROWING_TOPOLOGY_H

#include "instant_topology.h"

#include "island.h"

/// Base class for topologies than gan be "grown".
/**
 * Many topologies can be constructed incrementally, i.e. by adding one island after another.
 * Such topologies should inherit from this class.
 * Obviously, every "growing" topology is also an "instant creation" topology.
 */
class __PAGMO_VISIBLE growing_topology: public instant_topology {
	
	public:
		/// Virtual Destructor.
		virtual ~growing_topology() { };
		
		/// Add an island to the topology.
		/**
		 * \param[in] island The island to be added.
		 */
		virtual void push_back(const island &) = 0;
		
		/// \see instant_topology::isValidNumberOfNodes
		virtual bool isValidNumberOfNodes(const size_t& nodesCount) { return nodesCount >= 1; }
		
	protected:
		/// \see instant_topology::initialize
		virtual void do_initialize(const std::vector<island>& islands)
		{
			for(std::vector<island>::const_iterator i = islands.begin(); i != islands.end(); i++) {
				push_back(*i);
			}
		}
};

#endif
