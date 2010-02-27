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
#include <boost/python/operators.hpp>
#include <string>

#include "../../src/problem/base.h"
#include "../../src/problem/paraboloid.h"
#include "../../src/types.h"
#include "../exceptions.h"
#include "../utils.h"

using namespace boost::python;
using namespace pagmo;

template <class Problem>
static inline class_<Problem,bases<problem::base> > problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base> > retval(name,descr,init<const Problem &>());
	retval.def("__copy__", &Problem::clone);
	retval.def("__repr__", &Problem::human_readable);
	// Dimensions.
	retval.add_property("dimension", &Problem::get_dimension, "Global dimension.");
	retval.add_property("i_dimension", &Problem::get_i_dimension, "Integer dimension.");
	retval.add_property("c_dimension", &Problem::get_c_dimension, "Global constraints dimension.");
	retval.add_property("ic_dimension", &Problem::get_ic_dimension, "Inequality constraints dimension.");
	// Bounds.
	typedef void (Problem::*bounds_setter)(const decision_vector &);
	retval.add_property("lb",make_function(&Problem::get_lb,return_value_policy<copy_const_reference>()), bounds_setter(&Problem::set_lb), "Lower bounds.");
	retval.add_property("ub",make_function(&Problem::get_ub,return_value_policy<copy_const_reference>()), bounds_setter(&Problem::set_ub), "Upper bounds.");
	typedef void (Problem::*bounds_setter_value)(const double &, const double &);
	typedef void (Problem::*bounds_setter_vectors)(const decision_vector &, const decision_vector &);
	retval.def("set_bounds",bounds_setter_value(&Problem::set_bounds),"Set all bounds to the input values.");
	retval.def("set_bounds",bounds_setter_vectors(&Problem::set_bounds),"Set bounds to the input vectors.");
	retval.add_property("diameter",&Problem::get_diameter, "Problem's diameter.");
	// Useful operators.
	retval.def(self == self);
	retval.def(self != self);
	retval.def("is_compatible",&Problem::is_compatible,"Check compatibility with other problem.");
	// Comparisons.
	retval.def("compare_x",&Problem::compare_x,"Compare decision vectors.");
	retval.def("verify_x",&Problem::verify_x,"Check if decision vector is compatible with problem.");
	retval.def("compare_fc",&Problem::compare_fc,"Simultaneous fitness-constraint comparison.");
	// Constraints.
	typedef constraint_vector (Problem::*return_constraints)(const decision_vector &) const;
	retval.def("compare_constraints",&Problem::compare_constraints,"Compare constraint vectors.");
	retval.def("compute_constraints",return_constraints(&Problem::compute_constraints),"Compute and return constraint vector.");
	retval.def("feasibility_x",&Problem::feasibility_x,"Determine feasibility of decision vector.");
	retval.def("feasibility_c",&Problem::feasibility_c,"Determine feasibility of constraint vector.");
	// Fitness.
	typedef fitness_vector (Problem::*return_fitness)(const decision_vector &) const;
	retval.def("objfun",return_fitness(&Problem::objfun),"Compute and return fitness vector.");
	retval.def("compare_fitness",&Problem::compare_fitness,"Compare fitness vectors.");
	return retval;
}

BOOST_PYTHON_MODULE(_problem) {
	// Translate exceptions for this module.
	translate_exceptions();

	// Expose base problem class.
	class_<problem::base,boost::noncopyable>("__base", no_init);

	// Paraboloid problem.
	problem_wrapper<problem::paraboloid>("paraboloid","Multi-dimensional paraboloid miminisation.")
		.def(init<>())
		.def(init<const decision_vector &, const decision_vector &>());
}
