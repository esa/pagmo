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

#ifndef PAGMO_CHAIN_TOPOLOGY_H
#define PAGMO_CHAIN_TOPOLOGY_H

#include "graph_topology.h"

/// Chain topology (one-directional).
class __PAGMO_VISIBLE chain_topology: public graph_topology
{
	public:
		/// Constructor.
		chain_topology();
		/// Copy constructor.
		chain_topology(const chain_topology &);

		/// \see base_topology::clone
		virtual chain_topology *clone() const {
			return new chain_topology(*this);
		}

		/// \see base_topology::push_back
		virtual void push_back(const size_t&);

		/// \see base_topology::id_object()
		virtual std::string id_object() const {
			return id_name();
		}

	private:
		/// Tracks the identifier of the last inserted node.
		size_t	m_last;

		/// \see graph_topology::operator=
		chain_topology &operator=(const chain_topology &);
};

#endif
