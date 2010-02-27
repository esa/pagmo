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

#ifndef PAGMO_CARTWHEEL_TOPOLOGY_H
#define PAGMO_CARTWHEEL_TOPOLOGY_H

#include "../../../config.h"
#include "graph_topology.h"

/// Cartwheel topology - a two-way ring with every two opposite nodes connected.
class __PAGMO_VISIBLE cartwheel_topology: public graph_topology
{
	public:
		/// Constructor.
		cartwheel_topology();
		/// Copy constructor.
		cartwheel_topology(const cartwheel_topology &);

		/// \see base_topology::clone
		virtual cartwheel_topology *clone() const {
			return new cartwheel_topology(*this);
		}

		/// \see base_topology::push_back
		virtual void push_back(const size_t& id);

		/// \see base_topology::id_object()
		virtual std::string id_object() const {
			return id_name();
		}

	private:
		/// Tracks the id of the first inserted node on the ring.
		size_t	m_a;
		/// Tracks the id of the last inserted node on the ring.
		size_t	m_b;
		/// Tracks the id of the last node on the first half of the ring.
		size_t	m_c;
		/// Tracks the id of the first node on the second half of the ring.
		size_t	m_d;

		/// \see graph_topology::operator=
		cartwheel_topology &operator=(const cartwheel_topology &);
};

#endif
