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

// 12/01/2009: Initial version by Francesco Biscani.

#include "base_topology.h"

std::ostream &operator<<(std::ostream &os, const base_topology &topology)
{
	std::vector<size_t> graph_nodes = topology.get_nodes();

	os << topology.id_name() << std::endl;

	for (std::vector<size_t>::const_iterator it = graph_nodes.begin(); it != graph_nodes.end(); ++it) {
		//The node
		os << *it;

		//The edges
		std::vector<size_t> neighbours = topology.get_neighbours_out(*it);

		if (neighbours.size() > 0) {
			os << "->";
			for (size_t i = 0; i < neighbours.size(); ++i) {
				os << neighbours[i];
				if (i < neighbours.size() - 1) {
					os << ',';
				}
			}
		}
		os << std::endl;
	}
	return os;
}

