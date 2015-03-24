/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#include <string>

#include "../exceptions.h"
#include "base.h"
#include "rim.h"

namespace pagmo { namespace topology {

/// Default constructor.
rim::rim():base() {}

base_ptr rim::clone() const
{
	return base_ptr(new rim(*this));
}

void rim::connect(const vertices_size_type &)
{
	// Store frequently-used variables.
	const vertices_size_type t_size = get_number_of_vertices();
	pagmo_assert(t_size != 0);
	switch (t_size) {
	case 1: {
			// If the topology was empty, do not do anything.
			break;
		}
	case 2: {
			add_edge(0,1);
			add_edge(1,0);
			break;
		}
	case 3: {
			// Add edge to the center.
			add_edge(0,2);
			add_edge(2,0);
			// Add 1-2 connection.
			add_edge(1,2);
			add_edge(2,1);
			break;
		}
	case 4: {
			// Add edge to the center.
			add_edge(0,3);
			add_edge(3,0);
			// Add 1-3 and 3-2 connections.
			add_edge(1,3);
			add_edge(3,1);
			add_edge(2,3);
			add_edge(3,2);
			break;
		}
	default: {
			// Add edge to the center.
			add_edge(0,t_size - 1);
			add_edge(t_size - 1,0);
			// Remove connection (previous last)-first.
			remove_edge(t_size - 2,1);
			remove_edge(1,t_size - 2);
			// Add connection (previous last)-(new last).
			add_edge(t_size - 2,t_size - 1);
			add_edge(t_size - 1,t_size - 2);
			// Add connection (new last)-(first).
			add_edge(t_size - 1,1);
			add_edge(1,t_size - 1);
		}
	}
}

std::string rim::get_name() const
{
	return "Rim";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::topology::rim)
