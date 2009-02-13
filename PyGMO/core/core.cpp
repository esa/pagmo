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

#include "../../src/GOclasses/basic/archipelago.h"
#include "../../src/GOclasses/basic/ba_topology.h"
#include "../../src/GOclasses/basic/base_topology.h"
#include "../../src/GOclasses/basic/individual.h"
#include "../../src/GOclasses/basic/island.h"
#include "../../src/GOclasses/basic/no_topology.h"
#include "../../src/GOclasses/basic/one_way_ring_topology.h"
#include "../../src/GOclasses/basic/population.h"
#include "../../src/GOclasses/basic/ring_topology.h"
#include "../../src/GOclasses/problems/ClassicProblems.h"
#include "../../src/GOclasses/problems/TrajectoryProblems.h"
#include "../../src/GOclasses/problems/twodee_problem.h"
#include "../../src/exceptions.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;

template <class T>
static inline string Py_repr_vector(const vector<T> &v)
{
	ostringstream tmp;
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

template <class Container, class Item>
static inline void Py_set_item_from_ra(Container &c, int n, const Item &x)
{
	c[n] = x;
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

template <class Problem>
static inline class_<Problem,bases<GOProblem> > problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<GOProblem> > retval(name,descr,init<const Problem &>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__repr__", &Py_repr_from_stream<Problem>);
	retval.def("objfun", &Problem::objfun, "Objective function.");
	retval.def("set_lb", &Problem::set_lb, "Set lower bound.");
	retval.def("set_ub", &Problem::set_ub, "Set upper bound.");
	retval.add_property("dimension", &Problem::getDimension, "Dimension.");
	retval.add_property("id_name", &Problem::id_name, "Identification name.");
	return retval;
}

template <class Topology>
static inline class_<Topology,bases<base_topology> > topology_wrapper(const char *name, const char *descr)
{
	class_<Topology,bases<base_topology> > retval(name,descr,init<const Topology &>());
	retval.def("__copy__", &Py_copy_from_ctor<Topology>);
	retval.def("__repr__", &Py_repr_from_stream<Topology>);
	retval.add_property("id_name", &Topology::id_name, "Identification name.");
	return retval;
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(island_evolve_overloads, evolve, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(archipelago_evolve_overloads, evolve, 0, 1)

// Instantiate the core module.
BOOST_PYTHON_MODULE(_core)
{
	// Translate exceptions for this module.
	translate_exceptions();

	// Expose std::vector<double>.
	class_<vector<double> > class_vd("vector_double","std::vector<double>");
	class_vd.def("__repr__", &Py_repr_vector<double>);
	class_vd.def(vector_indexing_suite<vector<double> >());

	// Expose individual class.
	class_<Individual> class_ind("individual", "Individual.", init<const GOProblem &>());
	class_ind.def(init<const Individual &>());
	class_ind.def("__copy__", &Py_copy_from_ctor<Individual>);
	class_ind.def("__repr__", &Py_repr_from_stream<Individual>);
	class_ind.add_property("fitness", &Individual::getFitness, "Fitness.");
	class_ind.add_property("decision_vector", make_function(&Individual::getDecisionVector, return_value_policy<copy_const_reference>()),
		"Decision vector.");
	class_ind.add_property("velocity", make_function(&Individual::getVelocity, return_value_policy<copy_const_reference>()),
		"Velocity.");

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
	class_pop.def("__setitem__", &Py_set_item_from_ra<Population,Individual>);
	class_pop.def("__repr__", &Py_repr_from_stream<Population>);
	class_pop.add_property("problem", make_function(&problem_getter<GOProblem,Population>,return_value_policy<manage_new_object>()), "Problem.");
	class_pop.def("append", &Population::push_back, "Append individual at the end of the population.");
	class_pop.def("insert", &Population::insert, "Insert individual before index.");
	class_pop.def("mean", &Population::evaluateMean, "Evaluate mean.");
	class_pop.def("std", &Population::evaluateStd, "Evaluate std.");
	class_pop.def("best", &Population::extractBestIndividual, "Copy of best individual.");
	class_pop.def("worst", &Population::extractWorstIndividual, "Copy of worst individual.");

	// Expose base GOProblem class.
	class_<GOProblem, boost::noncopyable>("__go_problem", "Base GO problem", no_init);

	// Expose problem classes.
	// Trajectory problems.
	problem_wrapper<messengerfullProb>("messenger_full_problem", "Messenger full problem.").def(init<>());
	problem_wrapper<messengerProb>("messenger_problem", "Messenger problem.").def(init<>());
	problem_wrapper<cassini1Prob>("cassini1_problem", "Cassini1 problem.").def(init<>());
	// Classical problems.
	problem_wrapper<TestProb>("test_problem", "Test problem.").def(init<int>());
	problem_wrapper<rastriginProb>("rastrigin_problem", "Rastrigin problem.").def(init<int>());
	problem_wrapper<schwefelProb>("schwefel_problem", "Schwefel problem.").def(init<int>());
	problem_wrapper<ackleyProb>("ackley_problem", "Ackely problem.").def(init<int>());
	problem_wrapper<rosenbrockProb>("rosenbrock_problem", "Rosenbrock problem.").def(init<int>());
	problem_wrapper<lennardjonesProb>("lennardjones_problem", "Lennard-Jones problem.").def(init<int>());
	problem_wrapper<levyProb>("levy_problem", "Levy problem.").def(init<int>());
	// Twodee problem.
	problem_wrapper<twodee_problem>("twodee_problem", "Twodee problem.").def(init<int>()).def(init<int,const std::string &>());

	// Expose island.
	class_<island> class_island("island", "Island.", init<const GOProblem &, const go_algorithm &, int>());
	class_island.def(init<const GOProblem &, const go_algorithm &>());
	class_island.def(init<const island &>());
	class_island.def("__copy__", &Py_copy_from_ctor<island>);
	class_island.def("__delitem__", &island::erase);
	class_island.def("__getitem__", &island::operator[], "Get a copy of individual.");
	class_island.def("__len__", &island::size);
	class_island.def("__setitem__", &island::set);
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

	// Base topology.
	class_<base_topology, boost::noncopyable> class_bt("__base_topology", "Base topology.", no_init);

	// Topologies.
	topology_wrapper<no_topology>("no_topology", "No topology.").def(init<>());
	topology_wrapper<ring_topology>("ring_topology", "Ring topology.").def(init<const double &>());
	topology_wrapper<one_way_ring_topology>("one_way_ring_topology", "One way ring topology.").def(init<const double &>());
	topology_wrapper<ba_topology>("ba_topology", "BA model topology.").def(init<int, int, const double &>());

	// Expose archipelago.
	typedef island &(archipelago::*arch_get_island)(int);
	class_<archipelago> class_arch("archipelago", "Archipelago", init<const GOProblem &>());
	class_arch.def(init<const GOProblem &, const go_algorithm &, int, int>());
	class_arch.def(init<const GOProblem &, const base_topology &>());
	class_arch.def(init<const GOProblem &, const base_topology &, const go_algorithm &, int, int>());
	class_arch.def(init<const archipelago &>());
	class_arch.def("__copy__", &Py_copy_from_ctor<archipelago>);
	class_arch.def("__getitem__", arch_get_island(&archipelago::operator[]), return_internal_reference<>());
	class_arch.def("__len__", &archipelago::size);
	class_arch.def("__setitem__", &Py_set_item_from_ra<archipelago,island>);
	class_arch.def("__repr__", &Py_repr_from_stream<archipelago>);
	class_arch.add_property("topology", make_function(&topology_getter<base_topology,archipelago>, return_value_policy<manage_new_object>()),
		&archipelago::set_topology, "Topology.");
	class_arch.def("append", &archipelago::push_back, "Append island.");
	class_arch.add_property("problem", make_function(&problem_getter<GOProblem,archipelago>, return_value_policy<manage_new_object>()), "Problem.");
	class_arch.def("join", &archipelago::join, "Block until evolution on each island has terminated.");
	class_arch.add_property("busy", &archipelago::busy, "True if at least one island is evolving, false otherwise.");
	class_arch.def("evolve", &archipelago::evolve, archipelago_evolve_overloads());
	class_arch.def("evolve_t", &archipelago::evolve_t, "Evolve islands for an amount of time.");
}
