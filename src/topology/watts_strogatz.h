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

#ifndef PAGMO_TOPOLOGY_WATTS_STROGATZ_H
#define PAGMO_TOPOLOGY_WATTS_STROGATZ_H

#include <cstddef>
#include <string>

#include "../config.h"
#include "../rng.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Watts-Strogatz network model.
/**
 * \image html ws.png "Example of Watts-Strogatz model with N = 20, K = 4 and beta = 0.1."
 * \image latex ws.png "Example of Watts-Strogatz model with N = 20, K = 4 and beta = 0.1." width=8cm
 *
 * The Watts-Strogatz model is a ring lattice network in which forward edges are rewired with random probability. Such a network
 * has small-world properties, including short average path lengths and high clustering.
 *
 * In this implementation the graph grows dynamically by rewiring all the connections each time an island is added. Note that up to the the first K + 1
 * insertions, the topology will be fully connected. Afterwards, the topology will be a proper Watts-Strogatz model.
 *
 * Since the addition of a single element to the topology implies the rewiring of the whole topology, for archipelago objects of large size it is advisable
 * to build the topology outside the archipelago specifying the number of islands it will contain, and use archipelago::set_topology() to apply it to an existing
 * (and possibly unconnected) archipelago.
 *
 * @see http://en.wikipedia.org/wiki/Watts_and_Strogatz_model
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE watts_strogatz: public base
{
	public:
		watts_strogatz(int = 10, const double & = 0.05, int = 0);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void connect(const vertices_size_type &);
	private:
		void rewire();
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<std::size_t &>(m_k);
			ar & const_cast<double &>(m_beta);
			ar & m_drng;
			ar & m_urng;
		}
		const std::size_t	m_k;
		const double		m_beta;
		rng_double		m_drng;
		rng_uint32		m_urng;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::topology::watts_strogatz)

#endif
