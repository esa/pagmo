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

// 13/01/2009: Initial version by Marek Ruciński.

#include "torus_topology.h"

torus_topology::torus_topology():graph_topology(), m_in_first(0), m_in_last(0), m_out_first(0), m_out_last(0) {}

torus_topology::torus_topology(const torus_topology &r)
		:graph_topology(r), m_in_first(r.m_in_first), m_in_last(r.m_in_last), m_out_first(r.m_out_first), m_out_last(r.m_out_last) { }

torus_topology &torus_topology::operator=(const torus_topology &)
{
	pagmo_assert(false);
	return *this;
}

void torus_topology::push_back(const size_t& id)
{
	// Store frequently-used variables.
	const size_t t_size = get_number_of_nodes();

	//Add the new node to the graph
	add_node(id);

	if (t_size % 2) { // add to outer ring
		switch (t_size) {
		case 1: // the first node (in the outer ring)
			m_out_first = id;
			break;

		case 3: // the second node - connect to the first (in the outer ring)
			add_edge(id, m_out_last);
			add_edge(m_out_last, id);
			break;

		case 5: // the third node - connect to the first and last (in the outer ring)
			add_edge(id, m_out_last);
			add_edge(m_out_last, id);
			add_edge(id, m_out_first);
			add_edge(m_out_first, id);
			break;

		default: //all next - put between last and first (in the outer ring)
			remove_edge(m_out_last, m_out_first);
			remove_edge(m_out_first, m_out_last);
			add_edge(m_out_last, id);
			add_edge(id, m_out_last);
			add_edge(m_out_first, id);
			add_edge(id, m_out_first);
			break;
		}

		// Additionally, cross connect with the last node inserted into the inner ring.
		add_edge(m_in_last, id);
		add_edge(id, m_in_last);

		m_out_last = id;
	} else { //add to inner ring
		switch (t_size) {
		case 0: // the first node
			m_in_first = id;
			break;

		case 2: // the second node - connect to the first
			add_edge(id, m_in_last);
			add_edge(m_in_last, id);
			break;

		case 4: // the third node - connect to the first
			add_edge(id, m_in_last);
			add_edge(m_in_last, id);
			add_edge(id, m_in_first);
			add_edge(m_in_first, id);
			break;

		default: //all next - put between last and first
			remove_edge(m_in_last, m_in_first);
			remove_edge(m_in_first, m_in_last);
			add_edge(m_in_last, id);
			add_edge(id, m_in_last);
			add_edge(m_in_first, id);
			add_edge(id, m_in_first);
			break;
		}

		m_in_last = id;
	}
}
