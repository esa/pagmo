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

// 25/02/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_FULLY_CONNECTED_TOPOLOGY_H
#define PAGMO_FULLY_CONNECTED_TOPOLOGY_H

#include "graph_topology.h"

/// Fully-connected topology.
class __PAGMO_VISIBLE fully_connected_topology: public graph_topology
{
	public:
		/// Constructor.
		fully_connected_topology();

		/// Copy constructor
		fully_connected_topology(const fully_connected_topology &);

		/// \see base_topology::clone
		virtual fully_connected_topology *clone() const {
			return new fully_connected_topology(*this);
		}

		/// \see base_topology::push_back
		virtual void push_back(const size_t&);

		/// \see base_topology::id_object()
		virtual std::string id_object() const {
			return id_name();
		}

	private:
		/// \see graph_topology::operator=
		fully_connected_topology &operator=(const fully_connected_topology &);
};

#endif
