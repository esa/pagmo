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

#ifndef PAGMO_TOPOLOGY_WATTS_STROGATZ_H
#define PAGMO_TOPOLOGY_WATTS_STROGATZ_H

#include <cstddef>

#include "../config.h"
#include "../rng.h"
#include "base.h"

namespace pagmo { namespace topology {

/// Watts-Strogatz network model.
/**
 * \image html ws.png "Example of Watts-Strogatz model with N = 100, K = 10 and beta = 0.05."
 * \image latex ws_large.png "Example of Watts-Strogatz model with N = 100, K = 10 and beta = 0.05." width=6cm
 *
 * The Watts-Strogatz model is a ring lattice network in which forward edges are rewired with random probability. Such a network
 * has small-world properties, including short average path lengths and high clustering.
 *
 * In this implementation the graph grows dynamically by rewiring all the connections each time an island is added. Note that up to the the first K + 1
 * insertions, the topology will be fully connected. Afterwards, the topology will be a proper Watts-Strogatz model.
 *
 * @see http://en.wikipedia.org/wiki/Watts_and_Strogatz_model
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE watts_strogatz: public base
{
	public:
		watts_strogatz(int,const double &);
		base_ptr clone() const;
	protected:
		void connect(int);
	private:
		const std::size_t	m_k;
		const double		m_beta;
		rng_double		m_drng;
		rng_uint32		m_urng;

};

}}

#endif
