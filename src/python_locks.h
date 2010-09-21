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
 *   the Free Software Foundation; either version 3 of the License, or       *
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

#ifndef PAGMO_PYTHON_LOCKS_H
#define PAGMO_PYTHON_LOCKS_H

#include <Python.h>
#include <boost/thread/thread.hpp>
#include <boost/utility.hpp>

namespace pagmo {

class gil_state_lock: private boost::noncopyable
{
	public:
		gil_state_lock();
		~gil_state_lock();
	private:
		const boost::thread::id		m_cur_thread_id;
		PyGILState_STATE		m_gstate;
		static const boost::thread::id	m_main_thread_id;
};

class gil_releaser: private boost::noncopyable
{
	public:
		gil_releaser();
		~gil_releaser();
	private:
		PyThreadState *m_thread_state;
};

}

#endif
