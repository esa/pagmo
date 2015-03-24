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

#include <sstream>
#include <string>

#include "../exceptions.h"
#include "../rng.h"
#include "base.h"
#include "erdos_renyi.h"

namespace pagmo { namespace topology {

/// Constructor from probability.
/**
 * Construct an Erdős-Rényi graph topology with given probability parameter. Allowed values for the probability are in the [0,1] range.
 * Note that if the probability is null, the topology reduces to an unconnected topology, whereas if the probability is unitary
 * the topology reduces to a fully_connected topology.
 *
 * @param[in] prob probability parameter for the Erdős-Rényi model.
 */
erdos_renyi::erdos_renyi(const double &prob):base(),m_prob(prob),m_drng(rng_generator::get<rng_double>())
{
	if (prob < 0 || prob > 1) {
		pagmo_throw(value_error,"probability must be in the [0,1] range");
	}
}

base_ptr erdos_renyi::clone() const
{
	return base_ptr(new erdos_renyi(*this));
}

void erdos_renyi::connect(const vertices_size_type &n)
{
	for (std::pair<v_iterator,v_iterator> vertices = get_vertices(); vertices.first != vertices.second; ++vertices.first) {
		// Connect n bidirectionally to the other nodes with probability m_prob. Also, avoid to connect n with itself.
		if (*vertices.first != n && m_drng() < m_prob) {
			add_edge(n,*vertices.first);
			add_edge(*vertices.first,n);
		}
	}
}

std::string erdos_renyi::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "\tprobability: " << m_prob << '\n';
	return oss.str();
}

std::string erdos_renyi::get_name() const
{
	return "Erdos-Renyi";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::topology::erdos_renyi)
