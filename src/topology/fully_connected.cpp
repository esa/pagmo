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

#include <string>
#include <utility>

#include "fully_connected.h"

namespace pagmo { namespace topology {

/// Default constructor.
fully_connected::fully_connected():base() {}

base_ptr fully_connected::clone() const
{
	return base_ptr(new fully_connected(*this));
}

void fully_connected::connect(const vertices_size_type &n)
{
	for (std::pair<v_iterator,v_iterator> vertices = get_vertices(); vertices.first != vertices.second; ++vertices.first) {
		// Do not connect with self.
		if (n != *vertices.first) {
			add_edge(n,*vertices.first);
			add_edge(*vertices.first,n);
		}
	}
}

std::string fully_connected::get_name() const
{
	return "Fully connected";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::topology::fully_connected)
