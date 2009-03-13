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
#include <iostream>
#include <vector>

#include "../../Functions/rng/rng.h"
#include "individual.h"
#include "island.h"

/// Simple graph implementation to be used by topologies.
/**
 * \todo Rename this class.
 *
 * \todo I think this has to be re-designed. Interface of the ba_topology class should be extended with functions
 * for accessing vertices edges, and the graph_topology should be one of the available implementations of the underlying
 * graph. All other topologies (like ring, hypercube, etc.) are possible to be implemented with any underlying graph
 * implementation. At the moment I don't have an idea how to do this in an elegant way. Templates?
 */

class __PAGMO_VISIBLE graph_topology {
	
		///Stream output operator.
		friend std::ostream &operator<<(std::ostream &, const graph_topology &);
	
	protected:
		typedef boost::mutex mutex_type; ///< Mutex type abbreviation.
		typedef boost::lock_guard<mutex_type> lock_type; ///< Lock guard type abbreviation.
		
		typedef boost::unordered_map<size_t,std::vector<size_t> > tc_type; ///< Topology container type abbreviation.
		typedef tc_type::iterator tc_iterator; ///< Topology container iterator type abbreviation. \todo Maybe make the name a bit more verbose...
		
	public:
		///Default constructor.
		/** Creates a topology with no vertices nor edges */
		graph_topology();
		
		///Copy constructor.
		/** This strange thingy copies the structure but reinitialises the RNG... \todo Remove this ASAP */
		graph_topology(const graph_topology &);
	
	protected:
		mutable mutex_type	m_mutex; ///< Internal mutex to make certain operations thread-safe. \todo I'm not sure if it's still needed.
		tc_type				m_tc; ///< Graph structure - a map of lists of edges.
		rng_double			m_drng;	///< A RNG. \todo Move it elsewere.
	
	private:
		///Dummy assignment operator.
		/** Assignment is not a valid operation for topologies - throws exception when called. */
		graph_topology &operator=(const graph_topology &);
};

/// Stream output operator.
std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const graph_topology &);

#endif
