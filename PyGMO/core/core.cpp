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

// 27/12/2008: Initial version by Francesco Biscani.

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/def.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/module.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/utility.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/GOclasses/algorithms/go_algorithm.h"
#include "../../src/GOclasses/basic/archipelago.h"
#include "../../src/GOclasses/basic/base_topology.h"
#include "../../src/GOclasses/basic/individual.h"
#include "../../src/GOclasses/basic/island.h"
#include "../../src/GOclasses/basic/population.h"
#include "../../src/GOclasses/problems/GOproblem.h"
#include "../../src/GOclasses/basic/MigrationSelectionPolicy.h"
#include "../../src/GOclasses/basic/RandomMigrationSelectionPolicy.h"
#include "../../src/GOclasses/basic/MigrationReplacementPolicy.h"
#include "../../src/GOclasses/basic/RandomMigrationReplacementPolicy.h"
#include "../../src/exceptions.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;

template <class T>
static inline std::string Py_repr_vector(const std::vector<T> &v)
{
	std::ostringstream tmp;
	tmp << '<';
	const size_t size = v.size();
	for (size_t i = 0; i < size; ++i) {
		tmp << v[i];
		if (i + 1 != size) {
			tmp << ',';
		}
	}
	tmp << '>';
	return tmp.str();
}

template <class T, class C>
static inline T *problem_getter(const C &c)
{
	return c.problem().clone();
}

template <class T, class C>
static inline T *algorithm_getter(const C &c)
{
	return c.algorithm().clone();
}

template <class T, class C>
static inline T *topology_getter(const C &c)
{
	return c.topology().clone();
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(island_evolve_overloads, evolve, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(archipelago_evolve_overloads, evolve, 0, 1)

// Instantiate the core module.
BOOST_PYTHON_MODULE(_core)
{
	// Translate exceptions for this module.
	translate_exceptions();

	// Expose std::vector<double>.
	class_<std::vector<double> > class_vd("vector_double","std::vector<double>");
	class_vd.def("__repr__", &Py_repr_vector<double>);
	class_vd.def(vector_indexing_suite<std::vector<double> >());

	// Expose std::vector<int>.
	class_<std::vector<int> > class_vi("vector_int","std::vector<int>");
	class_vi.def("__repr__", &Py_repr_vector<int>);
	class_vi.def(vector_indexing_suite<std::vector<int> >());

	// Expose individual class.
	class_<Individual> class_ind("individual", "Individual.", init<const GOProblem &>());
	class_ind.def(init<const GOProblem &, const std::vector<double> &, const std::vector<double> &>());
	class_ind.def(init<const GOProblem &, const std::vector<double> &>());
	class_ind.def(init<const Individual &>());
	class_ind.def("__copy__", &Py_copy_from_ctor<Individual>);
	class_ind.def("__repr__", &Py_repr_from_stream<Individual>);
	class_ind.add_property("fitness", &Individual::getFitness, "Fitness.");
	class_ind.add_property("decision_vector", make_function(&Individual::getDecisionVector, return_value_policy<copy_const_reference>()),
		"Decision vector.");
	class_ind.add_property("velocity", make_function(&Individual::getVelocity, return_value_policy<copy_const_reference>()),
		"Velocity.");

	def("objfun_calls", &objfun_calls, "Number of times objective functions have been called.");

	// Expose population class.
	typedef const Individual &(Population::*pop_get_const)(int) const;
	class_<Population> class_pop("population", "Population.", init<const GOProblem &>());
	class_pop.def(init<const GOProblem &, int>());
	class_pop.def(init<const GOProblem &>());
	class_pop.def(init<const Population &>());
	class_pop.def("__copy__", &Py_copy_from_ctor<Population>);
	class_pop.def("__delitem__", &Population::erase);
	class_pop.def("__getitem__", pop_get_const(&Population::operator[]), return_value_policy<copy_const_reference>(), "Get a copy of individual.");
	class_pop.def("__len__", &Population::size);
	class_pop.def("__setitem__", &Population::setIndividual);
	class_pop.def("__repr__", &Py_repr_from_stream<Population>);
	class_pop.add_property("problem", make_function(&problem_getter<GOProblem,Population>,return_value_policy<manage_new_object>()), "Problem.");
	class_pop.def("append", &Population::push_back, "Append individual at the end of the population.");
	class_pop.def("insert", &Population::insert, "Insert individual before index.");
	class_pop.def("mean", &Population::evaluateMean, "Evaluate mean.");
	class_pop.def("std", &Population::evaluateStd, "Evaluate std.");
	class_pop.def("best", &Population::extractBestIndividual, return_value_policy<copy_const_reference>(), "Copy of best individual.");
	class_pop.def("worst", &Population::extractWorstIndividual, return_value_policy<copy_const_reference>(), "Copy of worst individual.");

	// Expose island.
	class_<island> class_island("island", "Island.", init<const GOProblem &, const go_algorithm &, int>());
	class_island.def(init<const GOProblem &, const go_algorithm &>());
	class_island.def(init<const island &>());
	class_island.def("__copy__", &Py_copy_from_ctor<island>);
	class_island.def("__delitem__", &island::erase);
	class_island.def("__getitem__", &island::operator[], "Get a copy of individual.");
	class_island.def("__len__", &island::size);
	class_island.def("__setitem__", &island::set_individual);
	class_island.def("__repr__", &Py_repr_from_stream<island>);
	class_island.def("append", &island::push_back, "Append individual at the end of the island.");
	class_island.def("insert", &island::insert, "Insert individual after index.");
	class_island.add_property("problem", make_function(&problem_getter<GOProblem,island>, return_value_policy<manage_new_object>()), "Problem.");
	class_island.add_property("algorithm", make_function(&algorithm_getter<go_algorithm,island>, return_value_policy<manage_new_object>()),
		&island::set_algorithm, "Algorithm.");
	class_island.add_property("population", &island::population, "Copy of population.");
	class_island.def("mean", &island::mean, "Evaluate mean.");
	class_island.def("std", &island::std, "Evaluate std.");
	class_island.def("best", &island::best, "Copy of best individual.");
	class_island.def("worst", &island::worst, "Copy of worst individual.");
	class_island.add_property("id", &island::id, "Identification number.");
	class_island.def("evolve", &island::evolve, island_evolve_overloads());
	class_island.def("evolve_t", &island::evolve_t, "Evolve for an amount of time.");
	class_island.def("join", &island::join, "Block until evolution has terminated.");
	class_island.add_property("busy", &island::busy, "True if island is evolving, false otherwise.");
	class_island.add_property("evo_time", &island::evo_time, "Total time spent evolving.");

	// Expose archipelago.
	typedef island &(archipelago::*arch_get_island)(int);
	class_<archipelago> class_arch("archipelago", "Archipelago", init<const GOProblem &>());
	class_arch.def(init<const GOProblem &, const go_algorithm &, int, int>());
	class_arch.def(init<const GOProblem &, const MigrationScheme&>());
	class_arch.def(init<const GOProblem &, const MigrationScheme&, const go_algorithm &, int, int, const MigrationSelectionPolicy&, const MigrationReplacementPolicy&>());
	class_arch.def(init<const archipelago &>());
	class_arch.def("__copy__", &Py_copy_from_ctor<archipelago>);
	class_arch.def("__getitem__", arch_get_island(&archipelago::operator[]), return_internal_reference<>());
	class_arch.def("__len__", &archipelago::size);
	class_arch.def("__setitem__", &archipelago::set_island);
	class_arch.def("__repr__", &Py_repr_from_stream<archipelago>);
	class_arch.add_property("migration_scheme", &archipelago::getMigrationScheme, &archipelago::setMigrationScheme, "Migration scheme.");
	class_arch.def("append", &archipelago::push_back, "Append island.");
	class_arch.add_property("problem", make_function(&problem_getter<GOProblem,archipelago>, return_value_policy<manage_new_object>()), "Problem.");
	class_arch.def("join", &archipelago::join, "Block until evolution on each island has terminated.");
	class_arch.add_property("busy", &archipelago::busy, "True if at least one island is evolving, false otherwise.");
	class_arch.def("evolve", &archipelago::evolve, archipelago_evolve_overloads());
	class_arch.def("evolve_t", &archipelago::evolve_t, "Evolve islands for an amount of time.");
	
	// Expose MigrationSelectionPolicy
	class_<MigrationSelectionPolicy, boost::noncopyable> class_MSP("__migration_selection_policy", "A migration selection policy.", no_init);
	
	// Expose RandomMigrationSelectionPolicy
	class_<RandomMigrationSelectionPolicy, bases<MigrationSelectionPolicy> > class_RMSP("random_migration_selection_policy", "A random migration selection policy.", init<optional<const uint32_t> >());	
	class_RMSP.def(init<const int&, optional<const uint32_t> >());
	class_RMSP.def(init<const int&, optional<const uint32_t> >());
	class_RMSP.def(init<const double&, optional<const uint32_t> >());
	
	
	// Expose MigrationReplacementPolicy
	class_<MigrationReplacementPolicy, boost::noncopyable> class_MRP("__migration_replacement_policy", "A migration replacement policy.", no_init);
	
	// Expose RandomMigrationReplacementPolicy
	class_<RandomMigrationReplacementPolicy, bases<MigrationReplacementPolicy> > class_RMRP("random_migration_replacement_policy", "A random migration replacement policy.", init<optional<const uint32_t> >());	
	class_RMSP.def(init<const int&, optional<const uint32_t> >());
	class_RMSP.def(init<const int&, optional<const uint32_t> >());
	class_RMSP.def(init<const double&, optional<const uint32_t> >());	
}
