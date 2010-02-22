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

#include <boost/numeric/conversion/cast.hpp>

#include "base.h"
#include "custom.h"

namespace pagmo { namespace topology {

/// Default constructor.
/**
 * Will call base::base().
 */
custom::custom():base() {}

/// Constructor from other topology.
/**
 * This constructor will copy the internal representation of any other topology (but of course not its push_back() logic).
 * Useful to take the snapshot of an existing topology and modify it manually for prototyping/experiments.
 *
 * @param[in] t topology to be copied.
 */
custom::custom(const base &t):base(t) {}

/// Check if a pair of island indices are adjacent.
/**
 * The direction of the edge must be n -> m. Will fail if either n or m is negative, or if either n or m is not in the topology.
 *
 * @param[in] n first island index.
 * @param[in] m second island index.
 *
 * @return true if the two islands are connected, false otherwise.
 */
bool custom::are_adjacent(int n, int m) const
{
	return base::are_adjacent(get_it(boost::numeric_cast<idx_type>(n)),get_it(boost::numeric_cast<idx_type>(m)));
}

/// Add an edge.
/**
 * Add an edge connecting index n to index m. Will fail if either n or m is negative, if either n or m is not in the topology or
 * if n and m are already connected.
 *
 * @param[in] n first index.
 * @param[in] m second index.
 */
void custom::add_edge(int n, int m)
{
	base::add_edge(get_it(boost::numeric_cast<idx_type>(n)),get_it(boost::numeric_cast<idx_type>(m)));
}

/// Remove an edge.
/**
 * Remove the edge connecting index n to index m. Will fail if either n or m is negative, if either n or m is not in the topology or
 * if n and m are not connected.
 *
 * @param[in] n first index.
 * @param[in] m second index.
 */
void custom::remove_edge(int n, int m)
{
	base::remove_edge(get_it(boost::numeric_cast<idx_type>(n)),get_it(boost::numeric_cast<idx_type>(m)));
}

/// Connect implementation.
/**
 * This class will not automatically add any connection during push_back operations.
 *
 * @param[in] n index to be connected.
 */
void custom::connect(int n)
{
	(void)n;
}

/// Remove all edges.
/**
 * Equivalent to base::remove_all_edges().
 */
void custom::remove_all_edges()
{
	base::remove_all_edges();
}

/// Add vertex.
/**
 * Equivalent to both base::add_vertex().
 *
 * @param[in] n index to be added to the topology.
 */
void custom::add_vertex(int n)
{
	base::add_vertex(n);
}

/// Number of adjacent vertices.
/**
 * Adjacent vertices are those connected from the interested vertex.
 *
 * @param[in] n interested index.
 *
 * @return number of adjacent vertices.
 */
custom::edges_size_type custom::get_num_adjacent_vertices(int n) const
{
	return base::get_num_adjacent_vertices(get_it(n));
}

/// Remove vertex.
/**
 * Equivalent to base::remove_vertex(). Will fail if n is negative or
 * if n is not present inside the topology.
 *
 * @param[in] n index to be removed from the topology.
 */
void custom::remove_vertex(int n)
{
	base::remove_vertex(get_it(n));
}

/// Clone method.
base_ptr custom::clone() const
{
	return base_ptr(new custom(*this));
}

}}
