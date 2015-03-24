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
#include "one_way_ring.h"

namespace pagmo { namespace topology {

/// Default constructor.
one_way_ring::one_way_ring():base(),m_first(0),m_last(0) {}

base_ptr one_way_ring::clone() const
{
	return base_ptr(new one_way_ring(*this));
}

void one_way_ring::connect(const vertices_size_type &n)
{
	// Store frequently-used variables.
	const vertices_size_type t_size = get_number_of_vertices();
	pagmo_assert(t_size != 0);
	switch (t_size) {
	case 1: {
			// If the topology was empty, just update the id of the first element.
			m_first = n;
			break;
		}
	case 2: {
			pagmo_assert(n != m_first);
			// Add connections to the only existing element.
			add_edge(m_first,n);
			add_edge(n,m_first);
			break;
		}
	default: {
			// The current last must be connected to the new one.
			remove_edge(m_last,m_first);
			add_edge(m_last,n);
			add_edge(n,m_first);
		}
	}
	// Update the id of the last island.
	m_last = n;
}

std::string one_way_ring::get_name() const
{
	return "One way ring";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::topology::one_way_ring)
