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

#ifndef PAGMO_LATTICE_TOPOLOGY_H
#define PAGMO_LATTICE_TOPOLOGY_H

#include "../../../config.h"
#include "graph_topology.h"

/// Lattice topology nodes are organised in a rectangle lattice with adjacent nodes connected.
/**
 * The topology has two parameters row_ratio and column_ratio, which describe the desired
 * aspect ratio of the lattice. The topology is constructed in such a way, that if the actual
 * row/columns ratio is smaller than or equal to the row_ratio/column ratio then a new row is added,
 * and when it's greater - a new column is added.
 */
class __PAGMO_VISIBLE lattice_topology: public graph_topology
{
	public:
		/// Constructor. Creates a square lattice
		lattice_topology();
		/// Constructor. Allows specification of the lattice aspect ratio.
		lattice_topology(size_t _row_ratio, size_t _column_ratio);
		/// Copy constructor.
		lattice_topology(const lattice_topology &);

		/// \see base_topology::clone
		virtual lattice_topology *clone() const {
			return new lattice_topology(*this);
		}

		/// \see base_topology::push_back
		virtual void push_back(const size_t& id);

		/// \see base_topology::id_object()
		virtual std::string id_object() const;

	private:
		const size_t row_ratio; ///< Numerator of the lattice aspect ratio.
		const size_t column_ratio; ///< Denominator of the lattice aspect ratio.

		std::vector<size_t> last_row; ///< Recently added row.
		std::vector<size_t> last_column; ///< Recently added column.
		size_t no_of_rows; ///< Number of <b>full</b> rows in the lattice.
		size_t no_of_columns; ///< Number of <b>full</b> columns in the lattice.
		size_t next_in_row; ///< Number of the nodes in the last (possibly incomplete) row.
		size_t next_in_column; ///< Number of the nodes in the last (possibly incomplete) column.

		/// \see graph_topology::operator=
		lattice_topology &operator=(const lattice_topology &);
};

#endif
