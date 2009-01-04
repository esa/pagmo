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
#include <vector>

#include "island.h"

class archipelago {
		typedef std::vector<island> container_type;
		typedef boost::mutex mutex_type;
		typedef boost::lock_guard<mutex_type> lock_type;
	public:
		typedef island value_type;
		archipelago() {};
		archipelago(const archipelago &);
		void push_back(const value_type &);
		size_t size() const;
		island get_island(int &) const;
	private:
		container_type		m_container;
		mutable mutex_type	m_mutex;
};

#endif
