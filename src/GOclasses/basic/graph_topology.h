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

// 22/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_GRAPH_TOPOLOGY_H
#define PAGMO_GRAPH_TOPOLOGY_H

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/unordered_map.hpp>
#include <vector>

#include "../../Functions/rng/rng.h"
#include "../../config.h"
#include "individual.h"
#include "island.h"

class __PAGMO_VISIBLE graph_topology {
	protected:
		typedef boost::mutex mutex_type;
		typedef boost::lock_guard<mutex_type> lock_type;
		// ic_type = individual container type.
		typedef boost::unordered_map<size_t,Individual> ic_type;
		typedef ic_type::iterator ic_iterator;
		// tc_type = topology container type.
		typedef boost::unordered_map<size_t,std::vector<size_t> > tc_type;
		typedef tc_type::iterator tc_iterator;
	public:
		graph_topology(const double &);
		graph_topology(const graph_topology &);
	protected:
		void pre_hook(island &);
		void post_hook(island &);
		mutable mutex_type	m_mutex;
		tc_type				m_tc;
		ic_type				m_ic;
		rng_double			m_drng;
		const double		m_prob;
	private:
		graph_topology &operator=(const graph_topology &);
};

#endif
