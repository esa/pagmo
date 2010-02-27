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

// 06/02/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_ONE_WAY_RING_TOPOLOGY_H
#define PAGMO_ONE_WAY_RING_TOPOLOGY_H

#include "../../../config.h"
#include "../island.h"
#include "graph_topology.h"

///Uni-directional ring topology.
class __PAGMO_VISIBLE one_way_ring_topology: public graph_topology
{
	public:
		/// Constructor.
		one_way_ring_topology();
		/// Copy constructor.
		one_way_ring_topology(const one_way_ring_topology &);

		/// \see base_topology::clone
		virtual one_way_ring_topology *clone() const {
			return new one_way_ring_topology(*this);
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
		/// Tracks the identifier of the last inserted node.
		size_t	m_last;

		/// \see graph_topology::operator=
		one_way_ring_topology &operator=(const one_way_ring_topology &);
};

#endif
