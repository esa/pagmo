/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

// Workaround for http://mail.python.org/pipermail/new-bugs-announce/2011-March/010395.html
#ifdef _WIN32
#include <cmath>
#endif
 
#include <Python.h>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>

#include "../../src/topologies.h"
#include "../utils.h"

using namespace boost::python;
using namespace pagmo;

template <class Topology>
static inline class_<Topology,bases<topology::base> > topology_wrapper(const char *name, const char *descr)
{
	class_<Topology,bases<topology::base> > retval(name,descr,init<const Topology &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Topology>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Topology>);
	retval.def_pickle(python_class_pickle_suite<Topology>());
	retval.def("cpp_loads", &py_cpp_loads<Topology>);
	retval.def("cpp_dumps", &py_cpp_dumps<Topology>);
	return retval;
}

// Wrappers for methods taking unsinged ints with safe conversion from int.
static inline bool topology_are_adjacent(const topology::base &t, int n, int m)
{
	return t.are_adjacent(boost::numeric_cast<topology::base::vertices_size_type>(n),boost::numeric_cast<topology::base::vertices_size_type>(m));
}

static inline bool topology_are_inv_adjacent(const topology::base &t, int n, int m)
{
	return t.are_inv_adjacent(boost::numeric_cast<topology::base::vertices_size_type>(n),boost::numeric_cast<topology::base::vertices_size_type>(m));
}

static inline std::vector<topology::base::vertices_size_type> topology_get_adjacent_vertices(const topology::base &t, int n)
{
	return t.get_v_adjacent_vertices(boost::numeric_cast<topology::base::vertices_size_type>(n));
}

static inline std::vector<topology::base::vertices_size_type> topology_get_inv_adjacent_vertices(const topology::base &t, int n)
{
	return t.get_v_inv_adjacent_vertices(boost::numeric_cast<topology::base::vertices_size_type>(n));
}

static inline topology::base::edges_size_type topology_get_num_adjacent_vertices(const topology::base &t, int n)
{
	return t.get_num_adjacent_vertices(boost::numeric_cast<topology::base::vertices_size_type>(n));
}

static inline topology::base::edges_size_type topology_get_num_inv_adjacent_vertices(const topology::base &t, int n)
{
	return t.get_num_inv_adjacent_vertices(boost::numeric_cast<topology::base::vertices_size_type>(n));
}

BOOST_PYTHON_MODULE(_topology) {
	common_module_init();
	typedef void (topology::base::*set_weight_edge)(const topology::base::vertices_size_type &, const topology::base::vertices_size_type&, double);
	typedef void (topology::base::*set_weight_vertex)(const topology::base::vertices_size_type &, double);
	typedef void (topology::base::*set_weight_all)(double);
	typedef double (topology::base::*get_weight_double)(const topology::base::vertices_size_type &, const topology::base::vertices_size_type&) const;

	class_<topology::base,boost::noncopyable>("_base",no_init)
		.def("__repr__", &topology::base::human_readable)
		.add_property("number_of_vertices",&topology::base::get_number_of_vertices)
		.add_property("number_of_edges",&topology::base::get_number_of_edges)
		.def("get_average_shortest_path_length",&topology::base::get_average_shortest_path_length,"Calculate average shortest path length.")
		.def("get_clustering_coefficient",&topology::base::get_clustering_coefficient,"Calculate the clustering coefficient.")
		.def("get_degree_distribution",&topology::base::get_degree_distribution,"Calculate the degree distribution.")
		.def("push_back",&topology::base::push_back,"Add vertex to the topology and connect it.")
		.def("set_weight",set_weight_edge(&topology::base::set_weight),"Set weight.")
		.def("set_weight",set_weight_vertex(&topology::base::set_weight),"Set weight.")
		.def("set_weight",set_weight_all(&topology::base::set_weight),"Set weight.")
		.def("get_weight",get_weight_double(&topology::base::get_weight),"Get weight.")
		.def("are_adjacent",&topology_are_adjacent,"Check whether two vertices are adjacent.")
		.def("are_inv_adjacent",&topology_are_inv_adjacent,"Check whether two vertices are inversely adjacent.")
		.def("get_adjacent_vertices",&topology_get_adjacent_vertices,"Return list of adjacent vertices.")
		.def("get_inv_adjacent_vertices",&topology_get_inv_adjacent_vertices,"Return list of inversely adjacent vertices.")
		.def("get_num_adjacent_vertices",&topology_get_num_adjacent_vertices,"Return number of adjacent vertices.")
		.def("get_num_inv_adjacent_vertices",&topology_get_num_inv_adjacent_vertices,"Return number of inversely adjacent vertices.");

	// Topologies.
	topology_wrapper<topology::barabasi_albert>("barabasi_albert", "Barabasi-Albert topology.").def(init<optional<int,int> >());
	topology_wrapper<topology::clustered_ba>("clustered_ba", "Clustered Barabasi-Albert topology.").def(init<optional<int,int,double> >());
	topology_wrapper<topology::ageing_clustered_ba>("ageing_clustered_ba", "Clustered Barabasi-Albert with Ageing topology.").def(init<optional<int,int,double,int> >());
	topology_wrapper<topology::custom>("custom", "Custom topology.")
		.def(init<const topology::base &>())
		.def("add_edge",&topology::custom::add_edge,"Add edge.")
		.def("remove_edge",&topology::custom::remove_edge,"Remove edge.")
		.def("remove_all_edges",&topology::custom::remove_all_edges,"Remove all edges.");

	topology_wrapper<topology::erdos_renyi>("erdos_renyi", "Erdos-Renyi topology.").def(init<optional<const double &> >());
	topology_wrapper<topology::fully_connected>("fully_connected", "Fully connected topology.");
	topology_wrapper<topology::ring>("ring", "Ring topology.");
	topology_wrapper<topology::hypercube>("hypercube", "Hypercube topology.");
	topology_wrapper<topology::one_way_ring>("one_way_ring", "One-way ring topology.");
	topology_wrapper<topology::unconnected>("unconnected", "Unconnected topology.");
	topology_wrapper<topology::watts_strogatz>("watts_strogatz", "Watts-Strogatz topology.").def(init<optional<int,const double &,int> >());
	topology_wrapper<topology::pan>("pan", "Pan graph topology.");
	topology_wrapper<topology::rim>("rim", "Wheel rim topology.");

	// Register to_python conversion from smart pointer.
	register_ptr_to_python<topology::base_ptr>();
}
