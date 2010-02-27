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
#include "one_way_ring.h"

namespace pagmo { namespace topology {

/// Default constructor.
one_way_ring::one_way_ring():base(),m_first(0),m_last(0) {}

/// Clone method.
base_ptr one_way_ring::clone() const
{
	return base_ptr(new one_way_ring(*this));
}

/// Connect method.
/**
 * Will insert the index into the ring topology, connecting it to the first index in the topology
 * and connecting the last index in the topology to it.
 */
void one_way_ring::connect(int n)
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
	default: {
			// The current last must be connected to the new one.
			const v_iterator it_first = get_it(m_first), it_n = get_it(n), it_last = get_it(m_last);
			remove_edge(it_last,it_first);
			add_edge(it_last,it_n);
			add_edge(it_n,it_first);
		}
	}
	// Update the id of the last island.
	m_last = n;
}

}}
