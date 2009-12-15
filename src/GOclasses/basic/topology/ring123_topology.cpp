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

#include "ring123_topology.h"

ring123_topology::ring123_topology():graph_topology(), a(0), b(0), c(0), d(0), e(0), f(0) {}

ring123_topology::ring123_topology(const ring123_topology &r)
		:graph_topology(r), a(r.a), b(r.b), c(r.c), d(r.d), e(r.e), f(r.f) {}

ring123_topology &ring123_topology::operator=(const ring123_topology &)
{
	pagmo_assert(false);
	return *this;
}

void ring123_topology::push_back(const size_t& id)
{
	// Store frequently-used variables.
	const size_t t_size = get_number_of_nodes();

	//Add the new node to the graph
	add_node(id);

	switch (t_size) {
	case 0:
		// If topology is empty, update the id of the first element.
		a = id;
		break;

	case 1:
		// Add connections to the only existing element.
		add_edge(id, a);
		add_edge(a, id);
		b = id;
		break;

	case 2:
		// Add new connections
		add_edge(a, id);
		add_edge(id, a);
		add_edge(b, id);
		add_edge(id, b);
		c = id;
		break;

	case 3:
		// Add new connections
		add_edge(a, id);
		add_edge(id, a);
		add_edge(b, id);
		add_edge(id, b);
		add_edge(c, id);
		add_edge(id, c);
		d = id;
		break;

	case 4:
		// Add new connections

		add_edge(a, id);
		add_edge(id, a);
		add_edge(b, id);
		add_edge(id, b);
		add_edge(c, id);
		add_edge(id, c);
		add_edge(d, id);
		add_edge(id, d);

		e = id;
		break;

	case 5:
		// Add new connections
		add_edge(a, id);
		add_edge(id, a);
		add_edge(b, id);
		add_edge(id, b);
		add_edge(c, id);
		add_edge(id, c);
		add_edge(d, id);
		add_edge(id, d);
		add_edge(e, id);
		add_edge(id, e);

		f = id;
		break;

	case 6:
		// Add new connections and update the nodes
		add_edge(a, id);
		add_edge(id, a);
		add_edge(b, id);
		add_edge(id, b);
		add_edge(c, id);
		add_edge(id, c);
		add_edge(d, id);
		add_edge(id, d);
		add_edge(e, id);
		add_edge(id, e);
		add_edge(f, id);
		add_edge(id, f);
		d = e;
		e = f;
		f = id;
		break;

	default:
		// Insert the new node between f and a.

		// First, put the node in between a and f. Former +1 link a-f becomes now a +2 link. Links b-f and a-e become +3.
		add_edge(f, id);
		add_edge(id, f);
		add_edge(a, id);
		add_edge(id, a);

		// Now drop connections b-e and d-a, and c-f as they are no longer +3 but +4, replacing them by edges c-id and b-id, e-id and and d-id
		remove_edge(b, e);
		remove_edge(e, b);
		remove_edge(c, f);
		remove_edge(f, c);
		remove_edge(a, d);
		remove_edge(d, a);
		add_edge(id, c);
		add_edge(c, id);
		add_edge(id, b);
		add_edge(b, id);
		add_edge(id, e);
		add_edge(e, id);
		add_edge(id, d);
		add_edge(d, id);

		// Finally, update the nodes
		d = e;
		e = f;
		f = id;
		break;
	}
}
