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

#ifndef PAGMO_HYPERCUBE_TOPOLOGY_H
#define PAGMO_HYPERCUBE_TOPOLOGY_H

#include "../../../config.h"
#include "graph_topology.h"

/// Hypercube topology (the number of dimensions grows automatically).
class __PAGMO_VISIBLE hypercube_topology: public graph_topology
{
	public:
		/// Constructor.
		hypercube_topology();
		/// Copy constructor.
		hypercube_topology(const hypercube_topology &);

		/// \see base_topology::clone
		virtual hypercube_topology *clone() const {
			return new hypercube_topology(*this);
		}

		/// \see base_topology::push_back
		virtual void push_back(const size_t& id);

		/// \see base_topology::id_object()
		virtual std::string id_object() const {
			return id_name();
		};

	private:
		std::vector<size_t> nodes; ///< All previously inserted nodes.

		/// \see graph_topology::operator=
		hypercube_topology &operator=(const hypercube_topology &);
};

#endif
