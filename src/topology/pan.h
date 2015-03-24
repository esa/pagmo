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

#ifndef PAGMO_TOPOLOGY_PAN_H
#define PAGMO_TOPOLOGY_PAN_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Pan graph topology.
/**
 * \image html pan_graph.png "Pan graph topologies."
 * \image latex pan_graph.png "Pan graph topologies." width=5cm
 *
 * The n-pan graph is the graph obtained by joining a cycle graph (i.e., a topology::ring topology) to a singleton graph with a bridge.
 * This implementation has bidirectional edges in the ring and a single edge towards the singleton, so that communication can flow from
 * the ring to the singleton, but not vice-versa. The index of the singleton vertex is always 0.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE pan: public base
{
	public:
		pan();
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
		}
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::topology::pan)

#endif
