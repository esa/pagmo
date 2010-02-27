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

#include "lattice_topology.h"
#include <sstream>

lattice_topology::lattice_topology()
		:graph_topology(), row_ratio(1), column_ratio(1), no_of_rows(0), no_of_columns(0), next_in_row(0), next_in_column(0) { }

lattice_topology::lattice_topology(size_t _row_ratio, size_t _column_ratio)
		:graph_topology(), row_ratio(_row_ratio), column_ratio(_column_ratio), no_of_rows(0), no_of_columns(0), next_in_row(0), next_in_column(0) { }


lattice_topology::lattice_topology(const lattice_topology &r)
		:graph_topology(r), row_ratio(r.row_ratio), column_ratio(r.column_ratio), last_row(r.last_row), last_column(r.last_column), no_of_rows(r.no_of_rows), no_of_columns(r.no_of_columns), next_in_row(r.next_in_row), next_in_column(r.next_in_column) { }

lattice_topology &lattice_topology::operator=(const lattice_topology &)
{
	pagmo_assert(false);
	return *this;
}

void lattice_topology::push_back(const size_t& id)
{
	// Store frequently-used variables.
	const size_t t_size = get_number_of_nodes();

	// Add node to the graph
	add_node(id);

	if (t_size == 0) { //the first node
		no_of_rows = no_of_columns = 1;

		last_row.push_back(id);
		last_column.push_back(id);
	} else { // all subsequent nodes
		if (no_of_rows * column_ratio <= no_of_columns * row_ratio) { // Add a new row
			//Connect the new node to the neighbour in the row "above"
			add_edge(last_row[next_in_row], id);
			add_edge(id, last_row[next_in_row]);

			if (next_in_row > 0) { // If no the first in row, connect also to the previous one in the same row
				add_edge(last_row[next_in_row - 1], id);
				add_edge(id, last_row[next_in_row - 1]);
			}

			//Store the node as the most recently added
			last_row[next_in_row] = id;
			//If the node is last in the row, it also belongs to the last column
			if (next_in_row + 1 == no_of_columns) {
				last_column.push_back(id);
			}
			//Increment the next node position
			++next_in_row;
			//If the row is complete, increase the number of full rows and wrap
			if (next_in_row == no_of_columns) {
				++no_of_rows;
				next_in_row = 0;
			}
		} else { // Add a new column
			//Connect the new node to the neighbour in the column "to the left"
			add_edge(last_column[next_in_column], id);
			add_edge(id, last_column[next_in_column]);

			if (next_in_column > 0) { // If not the first in column, connect also to the node above
				add_edge(last_column[next_in_column - 1], id);
				add_edge(id, last_column[next_in_column - 1]);
			}

			//Store the node as the most recently added
			last_column[next_in_column] = id;
			//If the node is last in the column, it also belongs to the last row
			if (next_in_column + 1 == no_of_rows) {
				last_row.push_back(id);
			}
			//Increment the next node position
			++next_in_column;
			//If the column is complete, increase the number of full columns and wrap
			if (next_in_column == no_of_rows) {
				++no_of_columns;
				next_in_column = 0;
			}
		}
	}
}

std::string lattice_topology::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << row_ratio << "_" << column_ratio;
	return tmp.str();
}
