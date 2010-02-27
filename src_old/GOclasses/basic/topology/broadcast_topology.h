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

// 06/02/2009: Initial version by Marek Ruciński.

#ifndef PAGMO_BROADCAST_TOPOLOGY_H
#define PAGMO_BROADCAST_TOPOLOGY_H

#include "../../../rng.h"
#include "graph_topology.h"

/// Broadcast topology (one node in the center, all others as leaves, bi-directional).
class __PAGMO_VISIBLE broadcast_topology: public graph_topology
{
	public:
		/// Constructor.
		broadcast_topology();
		/// Copy constructor.
		broadcast_topology(const broadcast_topology &);

		/// \see base_topology::clone
		virtual broadcast_topology *clone() const {
			return new broadcast_topology(*this);
		}

		/// \see base_topology::push_back
		virtual void push_back(const size_t&);

		/// \see base_topology::id_object()
		virtual std::string id_object() const {
			return id_name();
		}

	private:
		/// Tracks the identifier of the first inserted node.
		size_t	m_first;

		/// \see graph_topology::operator=
		broadcast_topology &operator=(const broadcast_topology &);
};

#endif
