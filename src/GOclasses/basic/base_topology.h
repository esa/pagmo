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
//#include "island.h"

class island;

/// Base class for topologies.
/**
 * This class defines the interface that all topologies must implement.
 * \todo Rename this class in order to keep to the standard conventions.
 */
class __PAGMO_VISIBLE base_topology {
	public:
		/// Constructor.
		base_topology() { };
		
		/// Virtual Destructor.
		virtual ~base_topology() { };
		
		/// Add an island to the topology.
		/**
		 * \param[in] island The island to be added.
		 */
		virtual void push_back(const island &) = 0;
		
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
