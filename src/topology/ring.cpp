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
#include "ring.h"

namespace pagmo { namespace topology {

/// Default constructor.
ring::ring():base(),m_first(0),m_last(0) {}

base_ptr ring::clone() const
{
	return base_ptr(new ring(*this));
}

void ring::connect(const vertices_size_type &n)
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
	case 3: {
			// Add new connections.
			add_edge(m_last,n);
			add_edge(n,m_last);
			add_edge(m_first,n);
			add_edge(n,m_first);
			break;
		}
	default: {
			// In general we must change the back connection of the first,
			// the forward connection of the current last, and add the new last
			// with proper connections.
			remove_edge(m_last,m_first);
			remove_edge(m_first,m_last);
			add_edge(m_last,n);
			add_edge(n,m_last);
			add_edge(m_first,n);
			add_edge(n,m_first);
		}
	}
	// Update the id of the last island.
	m_last = n;
}

std::string ring::get_name() const
{
	return "Ring";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::topology::ring)
