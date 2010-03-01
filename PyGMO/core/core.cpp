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
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <vector>

#include "../../src/algorithm/base.h"
#include "../../src/archipelago.h"
#include "../../src/island.h"
#include "../../src/migration/base_r_policy.h"
#include "../../src/migration/base_s_policy.h"
#include "../../src/population.h"
#include "../../src/problem/base.h"
#include "../../src/topology/base.h"
#include "../boost_python_container_conversions.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;
using namespace pagmo;

static inline problem::base_ptr problem_from_pop(const population &pop)
{
	return pop.problem().clone();
}

// Instantiate the core module.
BOOST_PYTHON_MODULE(_core)
{
	// Translate exceptions for this module.
	translate_exceptions();

	// Enable handy automatic conversions.
	to_tuple_mapping<std::vector<double> >();
	from_python_sequence<std::vector<double>,variable_capacity_policy>();
	to_tuple_mapping<std::vector<int> >();
	from_python_sequence<std::vector<int>,variable_capacity_policy>();
        to_tuple_mapping<std::vector<std::vector<double> > >();
        from_python_sequence<std::vector<std::vector<double> >,variable_capacity_policy>();
	to_tuple_mapping<std::vector<std::vector<int> > >();
	from_python_sequence<std::vector<std::vector<int> >,variable_capacity_policy>();

	// Expose population class.
	class_<population>("population", "Population class.", init<const problem::base &,optional<int> >())
		.def(init<const population &>())
		.def("__copy__", &Py_copy_from_ctor<population>)
		.def("__len__", &population::size)
		.def("__repr__", &population::human_readable)
		.def("problem",&problem_from_pop);

	// Expose island class.
	class_<island>("island", "Island class.", init<const problem::base &, const algorithm::base &,
		optional<int,double,const migration::base_s_policy &,const migration::base_r_policy &> >())
		.def(init<const island &>())
		.def("__copy__", &Py_copy_from_ctor<island>)
		.def("__len__", &island::get_size)
		.def("__repr__", &island::human_readable)
		.def("evolve", &island::evolve,"Evolve island n times.")
		.def("evolve_t", &island::evolve,"Evolve island for at least n milliseconds.")
		.def("join", &island::join,"Wait for evolution to complete.")
		.def("busy", &island::busy,"Check if island is evolving.");

	// Expose archipelago class.
	class_<archipelago>("archipelago", "Archipelago class.", init<const problem::base &, const algorithm::base &,
		int,int,optional<const topology::base &,archipelago::distribution_type,archipelago::migration_direction> >())
		.def(init<const archipelago &>())
		.def("__copy__", &Py_copy_from_ctor<archipelago>)
		.def("__len__", &archipelago::get_size)
		.def("__repr__", &archipelago::human_readable)
		.def("evolve", &archipelago::evolve,"Evolve archipelago n times.")
		.def("evolve_t", &archipelago::evolve_t,"Evolve archipelago for at least n milliseconds.")
		.def("join", &archipelago::join,"Wait for evolution to complete.")
		.def("busy", &archipelago::busy,"Check if archipelago is evolving.");
}
