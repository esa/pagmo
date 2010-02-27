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

#include "../exceptions.h"
#include "base.h"
#include "ring.h"

namespace pagmo { namespace topology {

/// Default constructor.
ring::ring():base(),m_first(0),m_last(0) {}

/// Clone method.
base_ptr ring::clone() const
{
	return base_ptr(new ring(*this));
}

/// Connect method.
/**
 * Will insert the index into the ring topology, connecting it to the first and the last indices in the topology
 * and connecting the last index in the topology to and from it.
 */
void ring::connect(int n)
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
			const v_iterator it_first = get_it(m_first), it_n = get_it(n);
			add_edge(it_first,it_n);
			add_edge(it_n,it_first);
			break;
		}
	case 3: {
			// Add new connections.
			const v_iterator it_first = get_it(m_first), it_n = get_it(n), it_last = get_it(m_last);
			add_edge(it_last,it_n);
			add_edge(it_n,it_last);
			add_edge(it_first,it_n);
			add_edge(it_n,it_first);
			break;
		}
	default: {
			const v_iterator it_first = get_it(m_first), it_n = get_it(n), it_last = get_it(m_last);
			// In general we must change the back connection of the first,
			// the forward connection of the current last, and add the new last
			// with proper connections.
			remove_edge(it_last,it_first);
			remove_edge(it_first,it_last);
			add_edge(it_last,it_n);
			add_edge(it_n,it_last);
			add_edge(it_first,it_n);
			add_edge(it_n,it_first);
		}
	}
	// Update the id of the last island.
	m_last = n;
}

}}

