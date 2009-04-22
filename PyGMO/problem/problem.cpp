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

// 13/02/2008: Initial version by Francesco Biscani.

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost_python_p_exceptions.h>
#include <string>

#include "../../src/GOclasses/problems/ClassicProblems.h"
#include "../../src/GOclasses/problems/GOproblem.h"
#include "../../src/GOclasses/problems/TrajectoryProblems.h"
#include "../../src/GOclasses/problems/twodee_problem.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;

template <class Problem>
static inline class_<Problem,bases<GOProblem> > problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<GOProblem> > retval(name,descr,init<const Problem &>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__repr__", &Py_repr_from_stream<Problem>);
	retval.def("objfun", &Problem::objfun, "Objective function.");
	retval.def("set_lb", &Problem::set_lb, "Set lower bound.");
	retval.def("set_ub", &Problem::set_ub, "Set upper bound.");
	retval.add_property("lb", make_function(&Problem::getLB,return_value_policy<copy_const_reference>()), &Problem::setLB, "Lower bounds.");
	retval.add_property("ub", make_function(&Problem::getUB,return_value_policy<copy_const_reference>()), &Problem::setUB, "Upper bounds.");
	retval.add_property("dimension", &Problem::getDimension, "Dimension.");
	retval.add_property("id_name", &Problem::id_name, "Identification name.");
	retval.add_property("id_object", &Problem::id_object, "Object identification name.");
	return retval;
}

BOOST_PYTHON_MODULE(_problem) {
	// Translate exceptions for this module.
	translate_p_exceptions();

	// Expose base GOProblem class.
	class_<GOProblem, boost::noncopyable>("__go_problem", "Base GO problem", no_init);

	// Expose problem classes.
	// Trajectory problems.
	problem_wrapper<messengerfullProb>("messenger_full", "Messenger full problem.").def(init<>());
	problem_wrapper<messengerProb>("messenger", "Messenger problem.").def(init<>());
	problem_wrapper<cassini1Prob>("cassini1", "Cassini1 problem.").def(init<>());
	problem_wrapper<laplaceProb>("laplace", "Laplace problem.").def(init<const std::vector<int> &>())
		.def("solution", &laplaceProb::solution, "Display solution relative to input decision vector.");
	// Classical problems.
	problem_wrapper<TestProb>("test", "Test problem.").def(init<int>());
	problem_wrapper<rastriginProb>("rastrigin", "Rastrigin problem.").def(init<int>());
	problem_wrapper<schwefelProb>("schwefel", "Schwefel problem.").def(init<int>());
	problem_wrapper<ackleyProb>("ackley", "Ackely problem.").def(init<int>());
	problem_wrapper<rosenbrockProb>("rosenbrock", "Rosenbrock problem.").def(init<int>());
	problem_wrapper<lennardjonesProb>("lennardjones", "Lennard-Jones problem.").def(init<int>());
	problem_wrapper<levyProb>("levy", "Levy problem.").def(init<int>());
	// Twodee problem.
	problem_wrapper<twodee_problem>("twodee", "Twodee problem.").def(init<int>()).def(init<int,const std::string &>());
}
