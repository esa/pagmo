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

#ifndef PAGMO_TOPOLOGY_BARABASI_ALBERT_H
#define PAGMO_TOPOLOGY_BARABASI_ALBERT_H

#include <cstddef>
#include <string>

#include "../config.h"
#include "../rng.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Barabási-Albert graph topology.
/**
 * \image html ba.png "Barabási-Albert network with m0 = 3, m = 2 and 100 vertices."
 * \image latex ba_large.png "Barabási-Albert network with m0 = 3, m = 2 and 100 vertices." width=7cm
 *
 * Topology based on the Barabási-Albert (BA) model for the generation of random undirected scale-free networks. The construction of this topology consists internally of
 * two phases:
 * - the first m0 elements added to the network constitute a kernel of nodes connected to each other with high probability;
 * - after the kernel is built, the next elements added to the network are connected randomly to m of the existing nodes; the probability
 *   of connection is biased linearly towards the most connected nodes.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 *
 * @see http://en.wikipedia.org/wiki/BA_model
 */
class __PAGMO_VISIBLE barabasi_albert: public base
{
	public:
		barabasi_albert(int m0 = 3, int m = 2);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void connect(const vertices_size_type &);
		std::string human_readable_extra() const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<std::size_t &>(m_m0);
			ar & const_cast<std::size_t &>(m_m);
			ar & m_drng;
			ar & m_urng;
		}  
		// Size of the kernel - the starting number of nodes.
		const std::size_t	m_m0;
		// Number of edges per newly-inserted node.
		const std::size_t	m_m;
		// Double random number generator
		rng_double		m_drng;
		// Integer random number generator.
		rng_uint32		m_urng;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::topology::barabasi_albert)

#endif
