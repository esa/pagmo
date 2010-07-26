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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <vector>

#include "../../src/algorithm/base.h"
#include "../../src/archipelago.h"
#include "../../src/island.h"
#include "../../src/keplerian_toolbox/keplerian_toolbox.h"
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

// Wrappers to make functions taking size_type as input take integers instead, with safety checks.
inline static island archipelago_get_island(const archipelago &a, int n)
{
	return a.get_island(boost::numeric_cast<archipelago::size_type>(n));
}

inline static void archipelago_set_algorithm(archipelago &archi, int n, const algorithm::base &a)
{
	archi.set_algorithm(boost::numeric_cast<archipelago::size_type>(n),a);
}

inline static population::individual_type population_get_individual(const population &pop, int n)
{
	return pop.get_individual(boost::numeric_cast<population::size_type>(n));
}

inline static void population_set_x(population &pop, int n, const decision_vector &x)
{
	pop.set_x(boost::numeric_cast<population::size_type>(n),x);
}

inline static void population_set_v(population &pop, int n, const decision_vector &v)
{
	pop.set_v(boost::numeric_cast<population::size_type>(n),v);
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
	to_tuple_mapping<std::vector<topology::base::vertices_size_type> >();
	from_python_sequence<std::vector<topology::base::vertices_size_type>,variable_capacity_policy>();
	to_tuple_mapping<std::vector<std::vector<double> > >();
	from_python_sequence<std::vector<std::vector<double> >,variable_capacity_policy>();
	to_tuple_mapping<std::vector<std::vector<int> > >();
	from_python_sequence<std::vector<std::vector<int> >,variable_capacity_policy>();
	to_tuple_mapping<std::vector<std::vector<topology::base::vertices_size_type> > >();
	from_python_sequence<std::vector<std::vector<topology::base::vertices_size_type> >,variable_capacity_policy>();
	// TODO: move this to kep toolbox exposition once it exists.
	to_tuple_mapping<kep_toolbox::array6D>();
	from_python_sequence<kep_toolbox::array6D,fixed_size_policy>();
	to_tuple_mapping<kep_toolbox::array3D>();
	from_python_sequence<kep_toolbox::array3D,fixed_size_policy>();
	enum_<kep_toolbox::epoch::type>("epoch_type")
		.value("MJD", kep_toolbox::epoch::MJD)
		.value("MJD2000", kep_toolbox::epoch::MJD2000)
		.value("JD", kep_toolbox::epoch::JD);
	class_<kep_toolbox::epoch>("epoch","Epoch class.",init<const double &,kep_toolbox::epoch::type>())
		.def(repr(self));
	class_<kep_toolbox::planet>("planet","Planet class.",init<const kep_toolbox::epoch&, const kep_toolbox::array6D&, const double &, const double &, const double &, const double &, optional<const std::string &> >())
		.def(repr(self));
	class_<kep_toolbox::planet_mpcorb,bases<kep_toolbox::planet> >("planet_mpcorb","Planet MPCORB class.",init<const std::string &>())
		.def("packed_date2epoch",&kep_toolbox::planet_mpcorb::packed_date2epoch)
		.staticmethod("packed_date2epoch");

	// Expose population class.
	class_<population>("population", "Population class.", init<const problem::base &,optional<int> >())
		.def(init<const population &>())
		.def("__copy__", &Py_copy_from_ctor<population>)
		.def("__getitem__", &population_get_individual)
		.def("__len__", &population::size)
		.def("__repr__", &population::human_readable)
		.add_property("problem",&problem_from_pop)
		.add_property("champion",make_function(&population::champion,return_value_policy<copy_const_reference>()))
		.def("get_best_idx",&population::get_best_idx,"Get index of best individual.")
		.def("get_worst_idx",&population::get_worst_idx,"Get index of worst individual.")
		.def("set_x", &population_set_x,"Set decision vector of individual at position n.")
		.def("set_v", &population_set_v,"Set velocity of individual at position n.")
		.def("push_back", &population::push_back,"Append individual with given decision vector at the end of the population.");

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
		optional<int,const double &,const migration::base_s_policy &,const migration::base_r_policy &> >())
		.def(init<const population &, const algorithm::base &,
		optional<const double &,const migration::base_s_policy &,const migration::base_r_policy &> >())
		.def(init<const island &>())
		.def("__copy__", &Py_copy_from_ctor<island>)
		.def("__len__", &island::get_size)
		.def("__repr__", &island::human_readable)
		.def("evolve", &island::evolve,"Evolve island n times.")
		.def("evolve_t", &island::evolve_t,"Evolve island for at least n milliseconds.")
		.def("join", &island::join,"Wait for evolution to complete.")
		.def("busy", &island::busy,"Check if island is evolving.")
		.def("is_blocking", &island::is_blocking,"Check if island is blocking.")
		.def("interrupt", &island::interrupt,"Interrupt evolution.")
		.add_property("problem",&island::get_problem)
		.add_property("algorithm",&island::get_algorithm,&island::set_algorithm)
		.add_property("population",&island::get_population);

	// Expose archipelago class.
	class_<archipelago>("archipelago", "Archipelago class.", init<const problem::base &, const algorithm::base &,
		int,int,optional<const topology::base &,archipelago::distribution_type,archipelago::migration_direction> >())
		.def(init<optional<archipelago::distribution_type,archipelago::migration_direction> >())
		.def(init<const topology::base &, optional<archipelago::distribution_type,archipelago::migration_direction> >())
		.def(init<const archipelago &>())
		.def("__copy__", &Py_copy_from_ctor<archipelago>)
		.def("__len__", &archipelago::get_size)
		.def("__repr__", &archipelago::human_readable)
		.def("__getitem__", &archipelago_get_island)
		.def("evolve", &archipelago::evolve,"Evolve archipelago n times.")
		.def("evolve_t", &archipelago::evolve_t,"Evolve archipelago for at least n milliseconds.")
		.def("join", &archipelago::join,"Wait for evolution to complete.")
		.def("interrupt", &archipelago::interrupt,"Interrupt evolution.")
		.def("busy", &archipelago::busy,"Check if archipelago is evolving.")
		.def("push_back", &archipelago::push_back,"Append island.")
		.def("set_algorithm", &archipelago_set_algorithm,"Set algorithm on island.")
		.def("is_blocking", &archipelago::is_blocking,"Check if archipelago is blocking.")
		.def("dump_migr_history", &archipelago::dump_migr_history)
		.def("clear_migr_history", &archipelago::clear_migr_history)
		.add_property("topology", &archipelago::get_topology, &archipelago::set_topology);

	// Archipelago's migration strategies.
	enum_<archipelago::distribution_type>("distribution_type")
		.value("point_to_point",archipelago::point_to_point)
		.value("broadcast",archipelago::broadcast);

	enum_<archipelago::migration_direction>("migration_direction")
		.value("source",archipelago::source)
		.value("destination",archipelago::destination);
}
