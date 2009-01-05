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

#ifndef PAGMO_ISLAND_H
#define PAGMO_ISLAND_H

#include <boost/scoped_ptr.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

#include "../../config.h"
#include "../algorithms/go_algorithm.h"
#include "../problems/GOproblem.h"
#include "population.h"

class archipelago;

class __PAGMO_VISIBLE island
{
		typedef boost::mutex mutex_type;
		typedef boost::lock_guard<mutex_type> lock_type;
		friend class archipelago;
	public:
		island(int, GOProblem &, const go_algorithm &);
		island(const island &);
		island &operator=(const island &);
		~island();
		Population get_pop() const;
		void evolve();
	private:
		void set_archipelago(archipelago *);
		struct evolver {
			evolver(island *i):m_i(i) {}
			void operator()();
			island *m_i;
		};
		Population				m_pop;
		boost::scoped_ptr<GOProblem>		m_gop;
		boost::scoped_ptr<go_algorithm>		m_goa;
		archipelago				*m_a;
		mutable mutex_type			m_mutex;
};

#endif
