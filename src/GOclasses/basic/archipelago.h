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

#ifndef PAGMO_ARCHIPELAGO_H
#define PAGMO_ARCHIPELAGO_H

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <list>

#include "../../config.h"
#include "../problems/GOproblem.h"
#include "../algorithms/go_algorithm.h"
#include "island.h"

class __PAGMO_VISIBLE archipelago {
		typedef boost::mutex mutex_type;
		typedef boost::lock_guard<mutex_type> lock_type;
		typedef std::list<island> container_type;
		friend std::ostream &operator<<(std::ostream &, const archipelago &);
	public:
		typedef container_type::iterator iterator;
		typedef container_type::const_iterator const_iterator;
		archipelago(const GOProblem &);
		archipelago(int, int, const GOProblem &, const go_algorithm &);
		archipelago(const archipelago &);
		void push_back(int, const go_algorithm &);
		const_iterator begin() const {return m_container.begin();}
		const_iterator end() const {return m_container.end();}
		iterator begin() {return m_container.begin();}
		iterator end() {return m_container.end();}
		size_t size() const;
		const GOProblem &problem() const;
		void join() const;
		bool busy() const;
		void evolve(int n = 1);
		void evolve_t(const double &);
	private:
		container_type						m_container;
		boost::scoped_ptr<const GOProblem>	m_gop;
		mutable mutex_type					m_mutex;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const archipelago &);

#endif
