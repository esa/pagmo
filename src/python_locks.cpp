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

#include <Python.h>
#include <boost/thread/thread.hpp>

#include "python_locks.h"

namespace pagmo {

const boost::thread::id gil_state_lock::m_main_thread_id = boost::this_thread::get_id();

// NOTE: this lock will check if it is being called by the main thread before doing anything.
gil_state_lock::gil_state_lock():m_cur_thread_id(boost::this_thread::get_id())
{
	if (m_cur_thread_id != m_main_thread_id) {
		m_gstate = PyGILState_Ensure();
	}
}

gil_state_lock::~gil_state_lock()
{
	if (m_cur_thread_id != m_main_thread_id) {
		PyGILState_Release(m_gstate);
	}
}

gil_releaser::gil_releaser()
{
	m_thread_state = PyEval_SaveThread();
}

gil_releaser::~gil_releaser()
{
	PyEval_RestoreThread(m_thread_state);
}

}
