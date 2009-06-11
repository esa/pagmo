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

// 13/02/2008: Initial version by Francesco Biscani.

#include <boost/cstdint.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>

#include "../../src/GOclasses/basic/topology/base_topology.h"
#include "../../src/GOclasses/basic/topology/graph_topology.h"
#include "../../src/GOclasses/basic/topology/ba_topology.h"
#include "../../src/GOclasses/basic/topology/fully_connected_topology.h"
#include "../../src/GOclasses/basic/topology/one_way_ring_topology.h"
#include "../../src/GOclasses/basic/topology/ring_topology.h"
#include "../../src/GOclasses/basic/topology/broadcast_topology.h"
#include "../../src/GOclasses/basic/topology/chain_topology.h"
#include "../../src/GOclasses/basic/topology/torus_topology.h"
#include "../../src/GOclasses/basic/topology/cartwheel_topology.h"
#include "../../src/GOclasses/basic/topology/lattice_topology.h"
#include "../../src/GOclasses/basic/topology/hypercube_topology.h"
#include "../../src/GOclasses/basic/topology/ring12_topology.h"
#include "../../src/GOclasses/basic/topology/ring123_topology.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;

template <class Topology>
static inline class_<Topology,bases<graph_topology> > topology_wrapper(const char *name, const char *descr)
{
	class_<Topology, bases<graph_topology> > retval(name, descr, init<const Topology &>());
	retval.def("__copy__", &Py_copy_from_ctor<Topology>);
	retval.def("__repr__", &Py_repr_from_stream<Topology>);
	retval.add_property("id_name", &Topology::id_name, "Identification name.");
	retval.add_property("id_object", &Topology::id_object, "Object identification name.");
	return retval;
}

BOOST_PYTHON_MODULE(_topology) {
	// Translate exceptions for this module.
	translate_exceptions();

	// Base topology.
	class_<base_topology, boost::noncopyable> class_bt("__base_topology", "Base topology.", no_init);
	
	// Graph topology implementation
	class_<graph_topology, bases<base_topology>, boost::noncopyable> graph_topo_class("__graph_topology", "Simple topology implementation.", no_init);

	// Topologies.
	topology_wrapper<ring_topology>("ring", "Ring topology.").def(init<>());
	topology_wrapper<one_way_ring_topology>("one_way_ring", "One way ring topology.").def(init<>());
	topology_wrapper<ba_topology>("ba", "BA model topology.").def(init<int, int, optional<boost::uint32_t> >());
	topology_wrapper<fully_connected_topology>("fully_connected", "Fully connected topology.").def(init<>());
	topology_wrapper<chain_topology>("chain", "Broadcast topology.").def(init<>());
	topology_wrapper<broadcast_topology>("broadcast", "Broadcast topology.").def(init<>());
	topology_wrapper<torus_topology>("torus", "Torus topology.").def(init<>());
	topology_wrapper<cartwheel_topology>("cartwheel", "Cartwheel topology.").def(init<>());
	topology_wrapper<lattice_topology>("lattice", "Lattice topology.").def(init<>()).def(init<size_t, size_t>());
	topology_wrapper<hypercube_topology>("hypercube", "Hypercube topology.").def(init<>());
	topology_wrapper<ring12_topology>("ring12", "+1+2 ring topology.").def(init<>());
	topology_wrapper<ring123_topology>("ring123", "+1+2+3 ring topology.").def(init<>());
}
