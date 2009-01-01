/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <boost/python/class.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/python/module.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/utility.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "../src/GOclasses/algorithms/ASA.h"
#include "../src/GOclasses/basic/individual.h"
#include "../src/GOclasses/basic/population.h"
#include "../src/GOclasses/problems/TrajectoryProblems.h"

using namespace boost::python;
using namespace std;

struct GOProblemWrap: GOProblem, wrapper<GOProblem>
{
	double objfun(const vector<double> &x)
	{
		return this->get_override("objfun")(x);
	}
};

struct go_algorithm_wrap: go_algorithm, wrapper<go_algorithm>
{
	Population evolve(const Population &pop, GOProblem &prob)
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

static inline void ie_translator(const index_error &ie) {
	PyErr_SetString(PyExc_IndexError, ie.what().c_str());
}

static inline void ve_translator(const value_error &ve) {
	PyErr_SetString(PyExc_ValueError, ve.what().c_str());
}

// Instantiate the PyGMO module.
BOOST_PYTHON_MODULE(_PyGMO)
{
	// Translate our C++ exceptions into Python exceptions.
	register_exception_translator<index_error>(ie_translator);
	register_exception_translator<value_error>(ve_translator);

	// Expose std::vector<double> and std::vector<size_t>.
	class_<vector<double> > class_vd("__base_vector_double","std::vector<double>");
	class_vd.def("__repr__", &Py_repr_vector<double>);
	class_vd.def(vector_indexing_suite<vector<double> >());

	class_<vector<size_t> > class_vs("__base_vector_size_t","std::vector<size_t>");
	class_vs.def("__repr__", &Py_repr_vector<size_t>);
	class_vs.def(vector_indexing_suite<vector<size_t> >());

	// Expose individual class.
	class_<Individual> class_ind("individual", "Individual.", init<GOProblem &>());
	class_ind.add_property("fitness", &Individual::getFitness, "Fitness.");
	class_ind.def("__repr__", &Py_repr_from_stream<Individual>);

	// Expose population class.
	class_<Population> class_pop("population", "Population.", init<>());
	class_pop.def(init<GOProblem &, int>());
	class_pop.def("__repr__", &Py_repr_from_stream<Population>);
	class_pop.def(vector_indexing_suite<Population>());
	class_pop.def("mean", &Population::evaluateMean, "Evaluate mean.");
	class_pop.def("std", &Population::evaluateStd, "Evaluate std.");
	class_pop.def("best", &Population::extractBestIndividual, return_internal_reference<>(), "Return best individual.");
	class_pop.def("worst", &Population::extractWorstIndividual, return_internal_reference<>(), "Return worst individual.");
	class_pop.def("extract_random_deme", &Population::extractRandomDeme, "Extract random deme.");

	// Expose base GOProblem class.
	class_<GOProblemWrap, boost::noncopyable> class_gop("go_problem", "Base GO problem", no_init);
	class_gop.def("objfun", pure_virtual(&GOProblem::objfun), "Objective function.");
	class_gop.add_property("dimension", &GOProblem::getDimension, "Dimension of the problem.");

	// Expose problem classes.
	class_<messengerfullProb, bases<GOProblem> > class_mfp("messenger_full_problem", "Messenger full problem.", init<>());

	// Expose base algorithm class.
	class_<go_algorithm_wrap, boost::noncopyable> class_goa("go_algorithm", "Base GO algorithm", no_init);
	class_goa.def("evolve", pure_virtual(&go_algorithm::evolve));

	// Expose algorithms.
	class_<ASAalgorithm, bases<go_algorithm> > class_asa("asa_algorithm", "ASA algorithm.", init<int, const GOProblem &, double, double>());
}
