/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

#include "../../config.h"
#include "base_topology.h"
#include "graph_topology.h"
#include "island.h"

/// Barabási-Albert scale-free topology.
/**
 * See http://en.wikipedia.org/wiki/BA_model and
 * http://www.nd.edu/~networks/Publication%20Categories/03%20Journal%20Articles/Physics/StatisticalMechanics_Rev%20of%20Modern%20Physics%2074,%2047%20(2002).pdf.
 */
class __PAGMO_VISIBLE ba_topology: public base_topology, public graph_topology {
	public:
		ba_topology(int, int, const double &);
		ba_topology(const ba_topology &);
		virtual ba_topology *clone() const {return new ba_topology(*this);}
		virtual void push_back(const island &);
		virtual void pre_evolution(island &);
		virtual void post_evolution(island &);
	private:
		ba_topology &operator=(const ba_topology &);
		// Starting number of nodes (m_0).
		const size_t	m_m_0;
		// Number of edges for newly-inserted nodes (m).
		const size_t	m_m;
};

#endif
