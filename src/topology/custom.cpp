/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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
#include <string>

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
 * Useful to take the snapshot of an existing topology and modify it manually.
 *
 * @param[in] t topology to be copied.
 */
custom::custom(const base &t):base(t) {}

/// Add an edge.
/**
 * Add an edge connecting index n to index m. Will fail if either n or m is negative, if either n or m is not in the topology or
 * if n and m are already connected.
 *
 * @param[in] n first index.
 * @param[in] m second index.
 * @param[in] migr_probability second index.
 */
void custom::add_edge(int n, int m, double migr_probability)
{
	base::add_edge(boost::numeric_cast<vertices_size_type>(n),boost::numeric_cast<vertices_size_type>(m));
	base::set_weight(n, m, migr_probability);
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
	base::remove_edge(boost::numeric_cast<vertices_size_type>(n),boost::numeric_cast<vertices_size_type>(m));
}

void custom::connect(const vertices_size_type &)
{}

/// Remove all edges.
/**
 * Equivalent to base::remove_all_edges().
 */
void custom::remove_all_edges()
{
	base::remove_all_edges();
}

base_ptr custom::clone() const
{
	return base_ptr(new custom(*this));
}

std::string custom::get_name() const
{
	return "Custom";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::topology::custom)
