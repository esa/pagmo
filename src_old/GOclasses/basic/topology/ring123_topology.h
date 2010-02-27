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

#ifndef PAGMO_RING123_TOPOLOGY_H
#define PAGMO_RING123_TOPOLOGY_H

#include "../../../config.h"
#include "graph_topology.h"

/// Bi-directional +1+2+3 ring topology
/**
 * In such a ring, every node is connected with a direct neigbour and his direct neighbour and his direct neighbour.
 */
class __PAGMO_VISIBLE ring123_topology: public graph_topology
{
	public:
		/// Constructor.
		ring123_topology();
		/// Copy constructor.
		ring123_topology(const ring123_topology &);

		/// \see base_topology::clone
		virtual ring123_topology *clone() const {
			return new ring123_topology(*this);
		}

		/// \see base_topology::push_back
		virtual void push_back(const size_t& id);

		/// \see base_topology::id_object()
		virtual std::string id_object() const {
			return id_name();
		}

	private:
		/// Tracks the id of the first tracked node.
		size_t	a;
		/// Tracks the id of the second tracked node.
		size_t	b;
		/// Tracks the id of the third tracked node.
		size_t	c;
		/// Tracks the id of the fourth tracked node.
		size_t	d;
		/// Tracks the id of the fifth tracked node.
		size_t	e;
		/// Tracks the id of the sixth tracked node.
		size_t	f;

		/// \see graph_topology::operator=
		ring123_topology &operator=(const ring123_topology &);
};

#endif
