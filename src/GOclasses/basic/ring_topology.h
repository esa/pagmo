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

// 13/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_RING_TOPOLOGY_H
#define PAGMO_RING_TOPOLOGY_H

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/unordered_map.hpp>
#include <vector>

#include "../../config.h"
#include "base_topology.h"
#include "graph_topology.h"
#include "individual.h"
#include "island.h"

class __PAGMO_VISIBLE ring_topology: public base_topology, public graph_topology {
	public:
		ring_topology(const double &);
		ring_topology(const ring_topology &);
		virtual ring_topology *clone() const {return new ring_topology(*this);}
		virtual void push_back(const island &);
		virtual void reset();
		virtual void pre_evolution(island &);
		virtual void post_evolution(island &);
	private:
		ring_topology &operator=(const ring_topology &);
		size_t	m_first;
		size_t	m_last;
};

#endif
