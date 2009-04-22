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

// 04/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_ARCHIPELAGO_H
#define PAGMO_ARCHIPELAGO_H

#include <list>

#include "../../config.h"
#include "../problems/GOproblem.h"
#include "../algorithms/go_algorithm.h"
#include "base_topology.h"
#include "island.h"
#include "py_container_utils.h"

/// Archipelago class.
class __PAGMO_VISIBLE archipelago: public py_container_utils<archipelago> {
		typedef std::list<island> container_type;
		typedef container_type::iterator iterator;
		typedef container_type::const_iterator const_iterator;
		friend std::ostream &operator<<(std::ostream &, const archipelago &);
		template <class T>
		friend class py_container_utils;
		friend class island;
// Work around behaviour of GCC < 4.1, which does not recognize
// friendship with classes defined inside friend classes.
#if GCC_VERSION < 401000
		friend class island::int_evolver;
		friend class island::t_evolver;
#endif
		const_iterator begin() const {return m_container.begin();}
		const_iterator end() const {return m_container.end();}
		iterator begin() {return m_container.begin();}
		iterator end() {return m_container.end();}
	public:
		// Ctors.
		archipelago(const GOProblem &);
		archipelago(const GOProblem &, const base_topology &);
		archipelago(const GOProblem &, const go_algorithm &, int, int);
		archipelago(const GOProblem &, const base_topology &, const go_algorithm &, int, int);
		archipelago(const archipelago &);
		const base_topology &topology() const;
		void set_topology(const base_topology &);
		const island &operator[](int) const;
		void set_island(int, const island &);
		void push_back(const island &);
		size_t size() const;
		const GOProblem &problem() const;
		void join() const;
		bool busy() const;
		void evolve(int n = 1);
		void evolve_t(const size_t &);
	private:
		archipelago &operator=(const archipelago &);
		void check_island(const island &) const;
		container_type						m_container;
		boost::scoped_ptr<const GOProblem>	m_gop;
		boost::scoped_ptr<base_topology>	m_top;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const archipelago &);

#endif
