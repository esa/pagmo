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

#ifndef PAGMO_TOPOLOGY_ONE_WAY_RING_H
#define PAGMO_TOPOLOGY_ONE_WAY_RING_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Uni-directional ring topology.
/**
 * \image html one_way_ring.png "Uni-directional ring topology example."
 * \image latex one_way_ring.png "Uni-directional ring topology example." width=4cm
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 */
class __PAGMO_VISIBLE one_way_ring: public base
{
	public:
		one_way_ring();
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void connect(const vertices_size_type &);
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_first;
			ar & m_last;
		}
		// Tracks the identifier of the first inserted vertex.
		vertices_size_type	m_first;
		// Tracks the identifier of the last inserted vertex.
		vertices_size_type	m_last;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::topology::one_way_ring)

#endif
