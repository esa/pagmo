/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

// 23/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_BA_TOPOLOGY_H
#define PAGMO_BA_TOPOLOGY_H

#include <boost/cstdint.hpp>

#include "../../../config.h"
#include "../../../Functions/rng/rng.h"
#include "graph_topology.h"

/// Barabasi-Albert scale-free network topology.
/**
 * See http://en.wikipedia.org/wiki/BA_model and
 * http://www.nd.edu/~networks/Publication%20Categories/03%20Journal%20Articles/Physics/StatisticalMechanics_Rev%20of%20Modern%20Physics%2074,%2047%20(2002).pdf.
 * \todo Rename this class.
 */

class __PAGMO_VISIBLE ba_topology: public graph_topology
{
	public:
		/// Constructor.
		/**
		 * Initialises the Barabasi-Albert topology generator.
		 * \param[in] m_0 Size of the kernel.
		 * \param[in] m Number of edges per new node.
		 * \param[in] optional random seed used to initialise the internal rng.
		 */
		ba_topology(int m_0, int m, boost::uint32_t seed = static_rng_uint32()());

		/// Copy constructor... \todo Change semantics and add a method to re-init an RNG?
		ba_topology(const ba_topology &);

		/// \see base_topology::clone
		virtual ba_topology *clone() const {
			return new ba_topology(*this);
		}

		/// \see base_topology::push_back
		virtual void push_back(const size_t& id);

		/// \see base_topology::id_object()
		virtual std::string id_object() const;

	private:
		/// \see graph_topology::operator=
		ba_topology &operator=(const ba_topology &);

		/// Size of the kernel - the starting number of nodes.
		const size_t			m_m_0;
		/// Number of edges per newly-inserted node.
		const size_t			m_m;

		/// Random number generator
		rng_double				drng;
		/// Seed with which rng was initialised.
		const boost::uint32_t	seed;
};

#endif
