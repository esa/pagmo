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
#include <boost/python/errors.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/python/module.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/utility.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/GOclasses/algorithms/ASA.h"
#include "../src/GOclasses/algorithms/DE.h"
#include "../src/GOclasses/algorithms/hs_algorithm.h"
#include "../src/GOclasses/basic/archipelago.h"
#include "../src/GOclasses/basic/individual.h"
#include "../src/GOclasses/basic/island.h"
#include "../src/GOclasses/basic/population.h"
#include "../src/GOclasses/problems/TrajectoryProblems.h"
#include "../src/exceptions.h"

using namespace boost::python;
using namespace std;

struct GOProblemWrap: GOProblem, wrapper<GOProblem>
{
#if defined ( __GNUC__ ) && GCC_VERSION < 400000
	GOProblemWrap(const size_t &s, const double *d1, const double *d2):GOProblem(s,d1,d2) {}
#endif
	double objfun(const vector<double> &x)
	{
		return this->get_override("objfun")(x);
	}
};

struct go_algorithm_wrap: go_algorithm, wrapper<go_algorithm>
{
	Population evolve(const Population &pop, const GOProblem &prob)
	{
		return this->get_override("evolve")(pop,prob);
	}
};

template <class T>
static inline string Py_repr_from_stream(const T &x)
{
	ostringstream tmp;
	tmp << x;
	return tmp.str();
}

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
static inline void Py_set_item_from_ra(Container &c, int n, const Item &x) {
	c[n] = x;
}

static inline void ie_translator(const index_error &ie) {
	PyErr_SetString(PyExc_IndexError, ie.what());
}

static inline void ve_translator(const value_error &ve) {
	PyErr_SetString(PyExc_ValueError, ve.what());
}

static inline void te_translator(const type_error &te) {
	PyErr_SetString(PyExc_TypeError, te.what());
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(island_evolve_overloads, evolve, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(archipelago_evolve_overloads, evolve, 0, 1)

// Instantiate the PyGMO module.
BOOST_PYTHON_MODULE(_PyGMO)
{
	// Translate our C++ exceptions into Python exceptions.
	register_exception_translator<index_error>(ie_translator);
	register_exception_translator<value_error>(ve_translator);
	register_exception_translator<type_error>(te_translator);

	// Expose std::vector<double> and std::vector<size_t>.
	class_<vector<double> > class_vd("__base_vector_double","std::vector<double>");
	class_vd.def("__repr__", &Py_repr_vector<double>);
	class_vd.def(vector_indexing_suite<vector<double> >());

	class_<vector<size_t> > class_vs("__base_vector_size_t","std::vector<size_t>");
	class_vs.def("__repr__", &Py_repr_vector<size_t>);
	class_vs.def(vector_indexing_suite<vector<size_t> >());

	// Expose individual class.
	class_<Individual> class_ind("individual", "Individual.", init<const GOProblem &>());
	class_ind.add_property("fitness", &Individual::getFitness, "Fitness.");
	class_ind.def("__repr__", &Py_repr_from_stream<Individual>);

	// Expose population class.
	typedef const Individual &(Population::*pop_get_const)(int) const;
	class_<Population> class_pop("population", "Population.", init<const GOProblem &>());
	class_pop.def(init<const GOProblem &, int>());
	class_pop.def(init<const GOProblem &>());
	class_pop.def("__delitem__", &Population::erase);
	class_pop.def("__getitem__", pop_get_const(&Population::operator[]), return_value_policy<copy_const_reference>(), "Get a copy of individual.");
	class_pop.def("__len__", &Population::size);
	class_pop.def("__setitem__", &Py_set_item_from_ra<Population,Individual>);
	class_pop.def("__repr__", &Py_repr_from_stream<Population>);
	class_pop.def("problem", &Population::problem, return_internal_reference<>(), "Return problem.");
	class_pop.def("append", &Population::push_back, "Append individual at the end of the population.");
	class_pop.def("insert", &Population::insert, "Insert individual after index.");
	class_pop.def("mean", &Population::evaluateMean, "Evaluate mean.");
	class_pop.def("std", &Population::evaluateStd, "Evaluate std.");
	class_pop.def("best", &Population::extractBestIndividual, "Copy of best individual.");
	class_pop.def("worst", &Population::extractWorstIndividual, "Copy of worst individual.");

	// Expose base GOProblem class.
	class_<GOProblemWrap, boost::noncopyable> class_gop("go_problem", "Base GO problem", no_init);
	class_gop.def("objfun", pure_virtual(&GOProblem::objfun), "Objective function.");
	class_gop.add_property("dimension", &GOProblem::getDimension, "Dimension of the problem.");
	class_gop.add_property("id_name", &GOProblem::id_name, "Identification name.");

	// Expose problem classes.
	class_<messengerfullProb, bases<GOProblem> > class_mfp("messenger_full_problem", "Messenger full problem.", init<>());
	class_<messengerProb, bases<GOProblem> > class_mp("messenger_problem", "Messenger problem.", init<>());

	// Expose base algorithm class.
	class_<go_algorithm_wrap, boost::noncopyable> class_goa("go_algorithm", "Base GO algorithm", no_init);
	class_goa.def("evolve", pure_virtual(&go_algorithm::evolve));
	class_goa.add_property("id_name", &go_algorithm::id_name, "Identification name.");

	// Expose algorithms.
	class_<ASAalgorithm, bases<go_algorithm> > class_asa("asa_algorithm", "ASA algorithm.", init<int, const double &, const double &>());
	class_<DEalgorithm, bases<go_algorithm> > class_de("de_algorithm", "DE algorithm.", init<int, const double &, const double &, int>());
	class_<hs_algorithm, bases<go_algorithm> > class_hs("hs_algorithm", "HS algorithm.", init<int, const double &, const double &, const double &>());

	// Expose island.
	class_<island> class_island("island", "Island.", init<const GOProblem &, const go_algorithm &, int>());
	class_island.def(init<const GOProblem &, const go_algorithm &>());
	class_island.def("__delitem__", &island::erase);
	class_island.def("__getitem__", &island::operator[], "Get a copy of individual.");
	class_island.def("__len__", &island::size);
	class_island.def("__setitem__", &island::set);
	class_island.def("__repr__", &Py_repr_from_stream<island>);
	class_island.def("append", &island::push_back, "Append individual at the end of the island.");
	class_island.def("insert", &island::insert, "Insert individual after index.");
	class_island.def("problem", &island::problem, return_internal_reference<>(), "Return problem.");
	class_island.def("algorithm", &island::algorithm, return_internal_reference<>(), "Return algorithm.");
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
	class_arch.def("__getitem__", arch_get_island(&archipelago::operator[]), return_internal_reference<>());
	class_arch.def("__len__", &archipelago::size);
	class_arch.def("__setitem__", &Py_set_item_from_ra<archipelago,island>);
	class_arch.def("__repr__", &Py_repr_from_stream<archipelago>);
	class_arch.def("append", &archipelago::push_back, "Append island.");
	class_arch.def("insert", &archipelago::insert, "Insert island after index.");
	class_arch.def("problem", &archipelago::problem, return_internal_reference<>(), "Return problem.");
	class_arch.def("join", &archipelago::join, "Block until evolution on each island has terminated.");
	class_arch.add_property("busy", &archipelago::busy, "True if at least one island is evolving, false otherwise.");
	class_arch.def("evolve", &archipelago::evolve, archipelago_evolve_overloads());
	class_arch.def("evolve_t", &archipelago::evolve_t, "Evolve islands for an amount of time.");
}
