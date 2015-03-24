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

#ifndef PAGMO_TOPOLOGY_CUSTOM_H
#define PAGMO_TOPOLOGY_CUSTOM_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Custom topology.
/**
 * This topology allows the user to manually build a topology by inserting nodes and creating connections between them. The connect() method
 * will leave new nodes unconnected. The intended use of this topology is to give the user the ability to create
 * quickly a custom topology without having to create another class, recompile, etc.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE custom: public base
{
	public:
		custom();
		custom(const base &);
		/** @name High-level graph manipulation for custom topologies. */
		//@{
		void add_edge(int,int, double = 1.0);
		void remove_edge(int,int);
		void remove_all_edges();
		//@}
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

BOOST_CLASS_EXPORT_KEY(pagmo::topology::custom)

#endif
