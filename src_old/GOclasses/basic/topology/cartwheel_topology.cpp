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

#include "cartwheel_topology.h"

cartwheel_topology::cartwheel_topology():graph_topology(), m_a(0), m_b(0), m_c(0), m_d(0) {}

cartwheel_topology::cartwheel_topology(const cartwheel_topology &r)
		:graph_topology(r), m_a(r.m_a), m_b(r.m_b), m_c(r.m_c), m_d(r.m_d) { }

cartwheel_topology &cartwheel_topology::operator=(const cartwheel_topology &)
{
	pagmo_assert(false);
	return *this;
}

void cartwheel_topology::push_back(const size_t& id)
{
	// Store frequently-used variables.
	const size_t t_size = get_number_of_nodes();

	//Add the new node to the graph
	add_node(id);

	switch (t_size) {
	case 0: // the first node
		m_a = id;
		break;

	case 1: // the second node - connect to the first
		add_edge(id, m_a);
		add_edge(m_a, id);
		m_c = id;
		break;

	case 2: // the third node - connect to the first and last (in the outer ring)
		add_edge(id, m_a);
		add_edge(m_a, id);
		add_edge(id, m_c);
		add_edge(m_c, id);
		m_d = id;
		break;

	case 3: // the fourth node - connect to the ring and add the diagonal connections
		add_edge(id, m_a);
		add_edge(m_a, id);
		add_edge(id, m_d);
		add_edge(m_d, id);

		// diagonal connections (note that m_a <-> m_d is already there)
		add_edge(m_c, id);
		add_edge(id, m_c);
		m_b = id;
		break;

	default: //all next
		if (t_size % 2) { //add between m_a and m_b
			remove_edge(m_c, m_d);
			remove_edge(m_d, m_c);

			add_edge(m_c, id);
			add_edge(id, m_c);
			add_edge(m_d, id);
			add_edge(id, m_d);

			//diagonal link
			add_edge(m_b, id);
			add_edge(id, m_b);

			// update m_c
			m_c = id;
		} else { //add between m_c and m_d and add the diagonal connection
			remove_edge(m_a, m_b);
			remove_edge(m_b, m_a);

			add_edge(m_a, id);
			add_edge(id, m_a);
			add_edge(m_b, id);
			add_edge(id, m_b);

			// update m_b
			m_b = id;
		}
		break;
	}
}
