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

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>

#include "../../src/topology/base.h"
#include "../../src/topology/ring.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;
using namespace pagmo;

template <class Topology>
static inline class_<Topology,bases<topology::base> > topology_wrapper(const char *name, const char *descr)
{
	class_<Topology,bases<topology::base> > retval(name,descr,init<const Topology &>());
	retval.def("__copy__", &Topology::clone);
	return retval;
}

BOOST_PYTHON_MODULE(_topology) {
	// Translate exceptions for this module.
	translate_exceptions();

	class_<topology::base,boost::noncopyable>("_base",no_init)
		.def("__repr__", &topology::base::human_readable);

	// Topologies.
	topology_wrapper<topology::ring>("ring", "Ring topology.").def(init<>());

	// Register to_python conversion from smart pointer.
	register_ptr_to_python<topology::base_ptr>();
}
