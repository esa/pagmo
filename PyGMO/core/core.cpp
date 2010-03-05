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

#define TRIVIAL_GETTER_SETTER(type1,type2,name) \
static inline type2 get_##name(const type1 &arg) \
{ \
	return arg.name; \
} \
static inline void set_##name(type1 &arg1, const type2 &arg2) \
{ \
	arg1.name = arg2; \
}

TRIVIAL_GETTER_SETTER(population::individual_type,decision_vector,cur_x);
TRIVIAL_GETTER_SETTER(population::individual_type,decision_vector,cur_v);
TRIVIAL_GETTER_SETTER(population::individual_type,fitness_vector,cur_f);
TRIVIAL_GETTER_SETTER(population::individual_type,constraint_vector,cur_c);
TRIVIAL_GETTER_SETTER(population::individual_type,decision_vector,best_x);
TRIVIAL_GETTER_SETTER(population::individual_type,fitness_vector,best_f);
TRIVIAL_GETTER_SETTER(population::individual_type,constraint_vector,best_c);

TRIVIAL_GETTER_SETTER(population::champion_type,decision_vector,x);
TRIVIAL_GETTER_SETTER(population::champion_type,fitness_vector,f);
TRIVIAL_GETTER_SETTER(population::champion_type,constraint_vector,c);

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
		.def("__getitem__", &population::get_individual,return_value_policy<copy_const_reference>())
		.def("__len__", &population::size)
		.def("__repr__", &population::human_readable)
		.add_property("problem",&problem_from_pop)
		.add_property("champion",make_function(&population::champion,return_value_policy<copy_const_reference>()))
		.def("set_x", &population::set_x,"Set decision vector of individual at position n.");

	// Individual and champion.
	class_<population::individual_type>("individual","Individual class.",init<>())
		.def("__repr__",&population::individual_type::human_readable)
		.add_property("cur_x",&get_cur_x,&set_cur_x)
		.add_property("cur_f",&get_cur_f,&set_cur_f)
		.add_property("cur_c",&get_cur_c,&set_cur_c)
		.add_property("cur_v",&get_cur_v,&set_cur_v)
		.add_property("best_x",&get_best_x,&set_best_x)
		.add_property("best_f",&get_best_f,&set_best_f)
		.add_property("best_c",&get_best_c,&set_best_c);

	class_<population::champion_type>("champion","Champion class.",init<>())
		.def("__repr__",&population::champion_type::human_readable)
		.add_property("x",&get_x,&set_x)
		.add_property("f",&get_f,&set_f)
		.add_property("c",&get_c,&set_c);

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
		.def("busy", &island::busy,"Check if island is evolving.")
		.def("interrupt", &island::interrupt,"Interrupt evolution.")
		.add_property("problem",&island::get_problem)
		.add_property("algorithm",&island::get_algorithm,&island::set_algorithm);

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
