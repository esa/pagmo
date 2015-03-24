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

#ifndef PAGMO_TOPOLOGY_UNCONNECTED_H
#define PAGMO_TOPOLOGY_UNCONNECTED_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Unconnected topology.
/**
 * \image html unconnected.png "Unconnected topology example."
 * \image latex unconnected.png "Unconnected topology example." width=4cm
 *
 * This topology will keep all islands disconnected from each other.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE unconnected: public base
{
	public:
		unconnected();
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

} }

BOOST_CLASS_EXPORT_KEY(pagmo::topology::unconnected)

 #endif
