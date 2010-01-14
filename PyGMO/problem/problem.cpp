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

// 13/02/2008: Initial version by Francesco Biscani.

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/def.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <string>

#include "../../src/GOclasses/problems/classic.h"
#include "../../src/GOclasses/problems/base.h"
#include "../../src/GOclasses/problems/trajectory.h"
#include "../../src/GOclasses/problems/earth_mars_lt.h"
#include "../../src/GOclasses/problems/earth_mars_lt2.h"
#include "../../src/GOclasses/problems/example.h"
#include "../../src/GOclasses/problems/twodee.h"
#include "../../src/GOclasses/problems/sp_test.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;
using namespace pagmo;

template <class Problem>
static inline class_<Problem,bases<problem::base> > problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base> > retval(name,descr,init<const Problem &>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__repr__", &Py_repr_from_stream<Problem>);
	retval.def("objfun", &Problem::objfun, "Objective function.");
	//retval.def("set_lb", &Problem::set_lb, "Set lower bound.");
	//retval.def("set_ub", &Problem::set_ub, "Set upper bound.");
	//retval.add_property("lb", make_function(&Problem::get_lb,return_value_policy<copy_const_reference>()), &Problem::set_lb, "Lower bounds.");
	//retval.add_property("ub", make_function(&Problem::get_ub,return_value_policy<copy_const_reference>()), &Problem::set_ub, "Upper bounds.");
	retval.add_property("dimension", &Problem::getDimension, "Dimension.");
	retval.add_property("id_name", &Problem::id_name, "Identification name.");
	retval.add_property("id_object", &Problem::id_object, "Object identification name.");
	return retval;
}

BOOST_PYTHON_MODULE(_problem) {
	// Translate exceptions for this module.
	translate_exceptions();

	// Expose base problem::base class.
	class_<problem::base, boost::noncopyable>("__go_problem", "Base GO problem", no_init);

	// Expose problem classes.
	// Trajectory problems from the GTOP database (http://www.esa.int/gsp/ACT/inf/op/globopt.htm)
	// MGA Problems
	problem_wrapper<problem::cassini1>("cassini1", "Cassini1 problem.").def(init<>());
	problem_wrapper<problem::gtoc1>("gtoc1", "Problem associated to the first edition of the GTOC competition.").def(init<>());
	//MGADSM Problems
	problem_wrapper<problem::messengerfull>("messenger_full", "Messenger full problem.").def(init<>());
	problem_wrapper<problem::messenger>("messenger", "Messenger problem.").def(init<>());
	problem_wrapper<problem::cassini2>("cassini2", "Cassini2 problem.").def(init<>());
	problem_wrapper<problem::rosetta>("rosetta", "Rosetta problem.").def(init<>());
	problem_wrapper<problem::sagas>("sagas", "Sagas problem.").def(init<>());
	problem_wrapper<problem::tandem>("tandem", "TandEM problem").def(init<const int>()).def(init<const int, const double>());

	//Miscellanea Trajectory Problems
	// MGADSM Problems
	problem_wrapper<problem::laplace>("laplace", "Laplace problem.").def(init<const std::vector<int> &>())
		.def("solution", &problem::laplace::solution, "Display solution relative to input decision vector.");
	//Low-Thrust Problems
	problem_wrapper<problem::earth_mars_lt>("earth_mars_lt", "Earth-Mars LT problem: impulsive transcription")
		.def(init<int, double, double, double>()).def("hr",&problem::earth_mars_lt::human_readable);
	problem_wrapper<problem::earth_mars_lt2>("earth_mars_lt2", "Earth-Mars LT problem: continuous transcription")
                .def(init<int, double, double, double>()).def("hr",&problem::earth_mars_lt2::human_readable)
                .def("visualize",&problem::earth_mars_lt2::visualize);
	
	// Classical problems.
	problem_wrapper<problem::test>("test1", "Example problem: minimization of y = sum x_i").def(init<int>());
	problem_wrapper<problem::example>("test2", "Example problem: minimization of y = x * x.").def(init<>());
	problem_wrapper<problem::rastrigin>("rastrigin", "Rastrigin problem.").def(init<int>());
	problem_wrapper<problem::schwefel>("schwefel", "Schwefel problem.").def(init<int>());
	problem_wrapper<problem::ackley>("ackley", "Ackely problem.").def(init<int>());
	problem_wrapper<problem::rosenbrock>("rosenbrock", "Rosenbrock problem.").def(init<int>());
	problem_wrapper<problem::lennardjones>("lennardjones", "Lennard-Jones problem.").def(init<int>());
	problem_wrapper<problem::levy>("levy", "Levy problem.").def(init<int>());
	// Stochastic Programming Problems (evolutionary robotics included)
	problem_wrapper<problem::twodee>("twodee", "Twodee problem.").def(init<int>()).def(init<int,const std::string &>());
	problem_wrapper<problem::inventory>("inventory", "Inventory problem.").def(init<int>());

	def("objfun_calls", &problem::objfun_calls, "Number of times objective functions have been called.");
	def("reset_objfun_calls", &problem::reset_objfun_calls, "Reset to zero the number of times objective functions have been called.");
}
