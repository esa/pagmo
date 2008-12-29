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

// Instantiate the PyGMO module.
BOOST_PYTHON_MODULE(_PyGMO)
{
	register_exception_translator<index_error>(ie_translator);

	class_<vector<double> > class_vd("__base_vector_double","std::vector<double>");
	class_vd.def("__repr__", &Py_repr_vector<double>);
	class_vd.def(vector_indexing_suite<vector<double> >());

	class_<Population> class_pop("population", "Population.", init<>());
	class_pop.def(init<GOProblem &, int>());
	class_pop.def("__repr__", &Py_repr_from_stream<Population>);
	class_pop.def(vector_indexing_suite<Population>());
	class_pop.add_property("mean", &Population::evaluateMean, "Evaluate mean.");
	class_pop.add_property("std", &Population::evaluateStd, "Evaluate std.");
	class_pop.def("best", &Population::extractBestIndividual, return_internal_reference<>(), "Return best individual.");
	class_pop.def("worst", &Population::extractWorstIndividual, return_internal_reference<>(), "Return worst individual.");

	class_<GOProblemWrap, boost::noncopyable> class_gop("goproblem", "Base GO problem", no_init);
	class_gop.def("objfun", pure_virtual(&GOProblem::objfun));
	class_gop.add_property("dimension", &GOProblem::getDimension, "Dimension of the problem.");

	class_<messengerfullProb, bases<GOProblem> > class_mfp("messenger_full_problem", "Messenger full problem.", init<>());

	class_<Individual> class_ind("individual", "Individual.", init<GOProblem &>());
	class_ind.add_property("fitness", &Individual::getFitness, "Fitness.");
	class_ind.def("__repr__", &Py_repr_from_stream<Individual>);
}
