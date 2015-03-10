/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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
#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/converter/registry.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/utility.hpp> // For boost::noncopyable.
#include <boost/array.hpp>
#include <sstream>
#include <vector>

#include "../../src/algorithm/base.h"
#include "../../src/archipelago.h"
#include "../../src/base_island.h"
#include "../../src/config.h"
#include "../../src/exceptions.h"
#include "../../src/migration/base_r_policy.h"
#include "../../src/migration/base_s_policy.h"
#include "../../src/migration/best_s_policy.h"
#include "../../src/migration/random_s_policy.h"
#include "../../src/migration/fair_r_policy.h"
#include "../../src/population.h"
#include "../../src/problem/base.h"
#include "../../src/topology/base.h"
#include "../boost_python_container_conversions.h"
#include "../utils.h"
#include "python_base_island.h"
#include "python_island.h"

#ifdef PAGMO_ENABLE_KEP_TOOLBOX
	#include "../../src/keplerian_toolbox/keplerian_toolbox.h"
#endif

using namespace boost::python;
using namespace pagmo;

static inline problem::base_ptr problem_from_pop(const population &pop)
{
	return pop.problem().clone();
}

static inline boost::python::tuple race_return_tuple(const population &pop, const population::size_type n_final,
									const unsigned int min_trials = 0,
									const unsigned int max_count = 1000,
									double delta = 0.05,
									const std::vector<population::size_type> &active_set = std::vector<population::size_type>(),
									const bool race_best = true,
									const bool screen_output = false) {
	std::pair<std::vector<pagmo::population::size_type>, unsigned int> res = pop.race(n_final,min_trials, max_count, delta, active_set, race_best ,screen_output);
	return boost::python::make_tuple(res.first,res.second);
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

TRIVIAL_GETTER_SETTER(population::individual_type,decision_vector,cur_x)
TRIVIAL_GETTER_SETTER(population::individual_type,decision_vector,cur_v)
TRIVIAL_GETTER_SETTER(population::individual_type,fitness_vector,cur_f)
TRIVIAL_GETTER_SETTER(population::individual_type,constraint_vector,cur_c)
TRIVIAL_GETTER_SETTER(population::individual_type,decision_vector,best_x)
TRIVIAL_GETTER_SETTER(population::individual_type,fitness_vector,best_f)
TRIVIAL_GETTER_SETTER(population::individual_type,constraint_vector,best_c)

TRIVIAL_GETTER_SETTER(population::champion_type,decision_vector,x)
TRIVIAL_GETTER_SETTER(population::champion_type,fitness_vector,f)
TRIVIAL_GETTER_SETTER(population::champion_type,constraint_vector,c)

// Wrappers to make functions taking size_type as input take integers instead, with safety checks.
inline static base_island_ptr archipelago_get_island(const archipelago &a, int n)
{
	return a.get_island(boost::numeric_cast<archipelago::size_type>(n));
}

inline static void archipelago_set_island(archipelago &a, int n, const base_island &isl)
{
	a.set_island(boost::numeric_cast<archipelago::size_type>(n),isl);
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

inline static void population_repair(population &pop, const int &idx, const algorithm::base_ptr &repair_algo)
{
	pop.repair(boost::numeric_cast<population::size_type>(idx),repair_algo);
}

struct __PAGMO_VISIBLE population_pickle_suite : boost::python::pickle_suite
{
	static boost::python::tuple getinitargs(const population &pop)
	{
		return boost::python::make_tuple(pop.problem().clone());
	}
	static boost::python::tuple getstate(const population &pop)
	{
		std::stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa << pop;
		return boost::python::make_tuple(ss.str(),pop.problem().clone());
	}
	static void setstate(population &pop, boost::python::tuple state)
	{
		if (len(state) != 2)
		{
			PyErr_SetObject(PyExc_ValueError,("expected 2-item tuple in call to __setstate__; got %s" % state).ptr());
			throw_error_already_set();
		}
		const std::string str = extract<std::string>(state[0]);
		std::stringstream ss(str);
		boost::archive::text_iarchive ia(ss);
		ia >> pop;
		const problem::base_ptr prob = boost::python::extract<problem::base_ptr>(state[1]);
		population_access::get_problem_ptr(pop) = prob->clone();
	}
};

template <class Island>
struct __PAGMO_VISIBLE island_pickle_suite : boost::python::pickle_suite
{
	static boost::python::tuple getinitargs(const Island &isl)
	{
		return boost::python::make_tuple(isl.get_algorithm(),isl.get_problem());
	}
	static boost::python::tuple getstate(const Island &isl)
	{
		std::stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa << isl;
		return boost::python::make_tuple(ss.str(),isl.get_algorithm(),isl.get_population());
	}
	static void setstate(Island &isl, boost::python::tuple state)
	{
		if (len(state) != 3)
		{
			PyErr_SetObject(PyExc_ValueError,("expected 3-item tuple in call to __setstate__; got %s" % state).ptr());
			throw_error_already_set();
		}
		const std::string str = extract<std::string>(state[0]);
		std::stringstream ss(str);
		boost::archive::text_iarchive ia(ss);
		ia >> isl;
		const algorithm::base_ptr algo = boost::python::extract<algorithm::base_ptr>(state[1]);
		isl.set_algorithm(*algo);
		const population pop = boost::python::extract<population>(state[2]);
		isl.set_population(pop);
	}
};

struct __PAGMO_VISIBLE archipelago_pickle_suite : boost::python::pickle_suite
{
	static boost::python::tuple getinitargs(const archipelago &)
	{
		return boost::python::make_tuple();
	}
	static boost::python::tuple getstate(const archipelago &archi)
	{
		std::stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa << archi;
		return boost::python::make_tuple(ss.str(),archi.get_islands());
	}
	static void setstate(archipelago &archi, boost::python::tuple state)
	{
		if (len(state) != 2)
		{
			PyErr_SetObject(PyExc_ValueError,("expected 2-item tuple in call to __setstate__; got %s" % state).ptr());
			throw_error_already_set();
		}
		const std::string str = extract<std::string>(state[0]);
		std::stringstream ss(str);
		boost::archive::text_iarchive ia(ss);
		ia >> archi;
		// Recover seaparately the islands.
		const std::vector<base_island_ptr> islands = extract<std::vector<base_island_ptr> >(state[1]);
		pagmo_assert(islands.size() == archi.get_size());
		for (std::vector<base_island_ptr>::size_type i = 0; i < islands.size(); ++i) {
			archi.set_island(i,*islands[i]);
		}
	}
};


#define REGISTER_CONVERTER(T,policy) \
{\
   boost::python::type_info info = boost::python::type_id<T >(); \
   const boost::python::converter::registration* reg = boost::python::converter::registry::query(info); \
   if (reg == NULL)  \
   {\
	to_tuple_mapping<T >();\
	from_python_sequence<T,policy>();\
   }\
   else if ((*reg).m_to_python == NULL)\
   {\
	to_tuple_mapping<T >();\
	from_python_sequence<T,policy>();\
   }\
}

// Instantiate the core module.
BOOST_PYTHON_MODULE(_core)
{

	common_module_init();
	typedef boost::array<double,2> array2D;
	//Register std converters to lists if not already registered by some other module
	REGISTER_CONVERTER(array2D,fixed_size_policy);
	REGISTER_CONVERTER(std::vector<double>, variable_capacity_policy);
	REGISTER_CONVERTER(std::vector<int>, variable_capacity_policy);
	REGISTER_CONVERTER(std::vector<topology::base::vertices_size_type>, variable_capacity_policy);
	REGISTER_CONVERTER(std::vector<std::vector<double> >, variable_capacity_policy);
	REGISTER_CONVERTER(std::vector<array2D>, variable_capacity_policy);
	REGISTER_CONVERTER(std::vector<std::vector<int> >, variable_capacity_policy);
	REGISTER_CONVERTER(std::vector<std::vector<topology::base::vertices_size_type> >, variable_capacity_policy);
	REGISTER_CONVERTER(std::vector<base_island_ptr>, variable_capacity_policy);
	REGISTER_CONVERTER(std::vector<pagmo::algorithm::base_ptr>, variable_capacity_policy);
	REGISTER_CONVERTER(std::vector<pagmo::problem::base_ptr>, variable_capacity_policy);
	
#ifdef PAGMO_ENABLE_KEP_TOOLBOX
	REGISTER_CONVERTER(std::vector<kep_toolbox::planet_ptr>, variable_capacity_policy);
#endif

	// Expose population class.

	typedef population::size_type (population::*get_best_1_idx)() const;
	typedef std::vector<population::size_type> (population::*get_best_N_idx)(const population::size_type& N) const;


	class_<population>("population", "Population class.", init<const problem::base &,optional<int, boost::uint32_t> >())
		.def(init<const population &>())
		.def("__copy__", &Py_copy_from_ctor<population>)
		.def("__deepcopy__", &Py_deepcopy_from_ctor<population>)
		.def("__getitem__", &population_get_individual)
		.def("__len__", &population::size)
		.def("__repr__", &population::human_readable)
		.add_property("problem",&problem_from_pop)
		.add_property("_problem_reference",make_function(&population::problem,return_value_policy<reference_existing_object>()))
		.add_property("champion",make_function(&population::champion,return_value_policy<copy_const_reference>()))
		.def("get_domination_list",&population::get_domination_list,return_value_policy<copy_const_reference>(), "Get the domination list for an indivdual")
		.def("compute_nadir",&population::compute_nadir, "Get the nadir objective vector")
		.def("compute_ideal",&population::compute_ideal, "Get the ideal objective vector")
		.def("compute_pareto_fronts",&population::compute_pareto_fronts, "Computes all Pareto fronts")
//		.def("get_crowding_d",&population::get_crowding_d, "returns crowding distance")
//		.def("update_pareto_information",&population::update_pareto_information, "updates crowding distance and front informations")
		.def("get_best_idx",get_best_1_idx(&population::get_best_idx),"Get index of best individual.")
		.def("get_best_idx",get_best_N_idx(&population::get_best_idx),"Get index of best N individual.")
		.def("get_worst_idx",&population::get_worst_idx,"Get index of worst individual.")
		.def("set_x", &population_set_x,"Set decision vector of individual at position n.")
		.def("set_v", &population_set_v,"Set velocity of individual at position n.")
		.def("push_back", &population::push_back,"Append individual with given decision vector at the end of the population.")
		.def("erase", &population::erase, "Erase individual at position")
		.def("mean_velocity", &population::mean_velocity, "Calculates the mean velocity across particles")
		.def("race", &race_return_tuple, "Race the individuals")
		.def("repair", &population_repair, "Repair the individual at the given index")
		.def("cpp_loads", &py_cpp_loads<population>)
		.def("cpp_dumps", &py_cpp_dumps<population>)
		.def_pickle(population_pickle_suite());

	// Individual and champion.
	class_<population::individual_type>("individual","Individual class.",init<>())
		.def("__repr__",&population::individual_type::human_readable)
		.def("__copy__", &Py_copy_from_ctor<population::individual_type>)
		.def("__deepcopy__", &Py_deepcopy_from_ctor<population::individual_type>)
		.add_property("cur_x",&get_cur_x,&set_cur_x)
		.add_property("cur_f",&get_cur_f,&set_cur_f)
		.add_property("cur_c",&get_cur_c,&set_cur_c)
		.add_property("cur_v",&get_cur_v,&set_cur_v)
		.add_property("best_x",&get_best_x,&set_best_x)
		.add_property("best_f",&get_best_f,&set_best_f)
		.add_property("best_c",&get_best_c,&set_best_c)
		.def("cpp_loads", &py_cpp_loads<population::individual_type>)
		.def("cpp_dumps", &py_cpp_dumps<population::individual_type>)
		.def_pickle(generic_pickle_suite<population::individual_type>());

	class_<population::champion_type>("champion","Champion class.",init<>())
		.def("__repr__",&population::champion_type::human_readable)
		.def("__copy__", &Py_copy_from_ctor<population::champion_type>)
		.def("__deepcopy__", &Py_deepcopy_from_ctor<population::champion_type>)
		.add_property("x",&get_x,&set_x)
		.add_property("f",&get_f,&set_f)
		.add_property("c",&get_c,&set_c)
		.def("cpp_loads", &py_cpp_loads<population::champion_type>)
		.def("cpp_dumps", &py_cpp_dumps<population::champion_type>)
		.def_pickle(generic_pickle_suite<population::champion_type>());

	// Base island class for Python implementation.
	class_<python_base_island, boost::noncopyable>("_base_island",init<const algorithm::base &, const problem::base &, optional<int,const migration::base_s_policy &,const migration::base_r_policy &> >())
		.def(init<const algorithm::base &, const population &, optional<const migration::base_s_policy &,const migration::base_r_policy &> >())
		.def("__repr__",&base_island::human_readable)
		.def("__len__", &base_island::get_size)
		.def("get_evolution_time", &base_island::get_evolution_time,"Gives the evolution time in milliseconds.")
		.def("evolve", &base_island::evolve,"Evolve island n times.")
		.def("evolve_t", &base_island::evolve_t,"Evolve island for at least n milliseconds.")
		.def("join", &base_island::join,"Wait for evolution to complete.")
		.def("busy", &base_island::busy,"Check if island is evolving.")
		.def("interrupt", &base_island::interrupt,"Interrupt evolution.")
		.def("set_x", &base_island::set_x, "Assigns a decision vector to the i-th individual of the island population")
		.def("set_v", &base_island::set_x, "Assigns a velocity vector to the i-th individual of the island population")
		.add_property("problem",&base_island::get_problem)
		.add_property("algorithm",&base_island::get_algorithm,&island::set_algorithm)
		.add_property("population",&base_island::get_population, &base_island::set_population)
		.add_property("s_policy",&base_island::get_s_policy)
		.add_property("r_policy",&base_island::get_r_policy)
		// Virtual methods.
		.def("get_name", &base_island::get_name,&python_base_island::default_get_name)
		.def("_perform_evolution",&python_base_island::py_perform_evolution)
		.def_pickle(python_base_island_pickle_suite());

	// Local island class.
	class_<python_island,bases<base_island> >("local_island", "Local island class.",init<const algorithm::base &, const problem::base &, optional<int,const migration::base_s_policy &,const migration::base_r_policy &> >())
		.def(init<const algorithm::base &, const population &, optional<const migration::base_s_policy &,const migration::base_r_policy &> >())
		.def(init<const python_island &>())
		.def("__copy__", &Py_copy_from_ctor<python_island>)
		.def("__deepcopy__", &Py_deepcopy_from_ctor<python_island>)
		.def("cpp_loads", &py_cpp_loads<python_island>)
		.def("cpp_dumps", &py_cpp_dumps<python_island>)
		.def("is_pythonic", &python_island::is_pythonic)
		.def_pickle(island_pickle_suite<python_island>());

	// Register to_python conversion from smart pointer.
	register_ptr_to_python<base_island_ptr>();

	// Expose archipelago class.
	class_<archipelago>("archipelago", "Archipelago class.", init<const algorithm::base &, const problem::base &,
		int,int,optional<const topology::base &,archipelago::distribution_type,archipelago::migration_direction> >())
		.def(init<optional<archipelago::distribution_type,archipelago::migration_direction> >())
		.def(init<const topology::base &, optional<archipelago::distribution_type,archipelago::migration_direction> >())
		.def(init<archipelago::distribution_type, archipelago::migration_direction>())
		.def(init<const archipelago &>())
		.def("__copy__", &Py_copy_from_ctor<archipelago>)
		.def("__deepcopy__", &Py_deepcopy_from_ctor<archipelago>)
		.def("__len__", &archipelago::get_size)
		.def("__repr__", &archipelago::human_readable)
		.def("__getitem__", &archipelago_get_island)
		.def("__setitem__", &archipelago_set_island)
		.def("get_islands", &archipelago::get_islands)
		.def("set_seeds", &archipelago::set_seeds)
		.def("evolve", &archipelago::evolve,"Evolve archipelago *n* times.",boost::python::args("n"))
		.def("evolve_batch", &archipelago::evolve_batch,"Evolve archipelago *n* times in batches of *b* islands.",boost::python::args("n","b"))
		.def("evolve_t", &archipelago::evolve_t,"Evolve archipelago for at least *n* milliseconds.",boost::python::args("n"))
		.def("join", &archipelago::join,"Wait for evolution to complete.")
		.def("interrupt", &archipelago::interrupt,"Interrupt evolution.")
		.def("busy", &archipelago::busy,"Check if archipelago is evolving.")
		.def("push_back", &archipelago::push_back,"Append island.")
		.def("set_algorithm", &archipelago_set_algorithm,"Set algorithm on island.")
		.def("dump_migr_history", &archipelago::dump_migr_history)
		.def("clear_migr_history", &archipelago::clear_migr_history)
		.def("cpp_loads", &py_cpp_loads<archipelago>,
			"Load C++ serialized representation from string *str*.\n\n"
			":Parameters:\n"
			"  str: string\n",
			boost::python::args("str"))
		.def("cpp_dumps", &py_cpp_dumps<archipelago>,
			"Dump C++ serialized representation.\n\n"
			":Returns:\n"
			"   string representing the serialized C++ representation\n"
		)
		.add_property("topology", &archipelago::get_topology, &archipelago::set_topology,"Topology property.")
		.add_property("distribution_type", &archipelago::get_distribution_type, &archipelago::set_distribution_type, "Distribution type property.")
		.def_pickle(archipelago_pickle_suite());

	// Archipelago's migration strategies.
	enum_<archipelago::distribution_type>("distribution_type")
		.value("point_to_point",archipelago::point_to_point)
		.value("broadcast",archipelago::broadcast);

	enum_<archipelago::migration_direction>("migration_direction")
		.value("source",archipelago::source)
		.value("destination",archipelago::destination);
}
