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

#ifndef PAGMO_TOPOLOGY_ERDOS_RENYI_H
#define PAGMO_TOPOLOGY_ERDOS_RENYI_H


#include <string>

#include "../config.h"
#include "../rng.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Erdős-Rényi graph topology.
/**
 * \image html erdos_renyi.png "Erdős-Rényi G(n,p) model with n = 100, p = 0.02."
 * \image latex erdos_renyi_large.png "Erdős-Rényi G(n,p) model with n = 100, p = 0.02." width=7cm
 *
 * In the Erdős-Rényi \f$ G(n,p) \f$ model (ER), a graph with \f$ n \f$ vertices is constructed by connecting vertices randomly,
 * so that each possible edge is included with probability \f$ p \f$ (independent from the presence or absence of any other edge
 * in the graph). The expected number of edges in \f$ G(n,p) \f$ is \f$ {n \choose 2} p \f$.
 *
 * In this implementation, each time an island is added to the topology each new possible bidirectional edge to and from the new island
 * is created with probability \f$ p \f$.
 *
 * @see http://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE erdos_renyi: public base
{
	public:
		erdos_renyi(const double &prob = 0.01);
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
			ar & const_cast<double &>(m_prob);
			ar & m_drng;
		}  
		const double	m_prob;
		rng_double	m_drng;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::topology::erdos_renyi)

#endif
