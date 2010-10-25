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
#include <boost/python/enum.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>
#include <cstddef>
#include <string>

#include "../../src/exceptions.h"
#include "../../src/keplerian_toolbox/keplerian_toolbox.h"
#include "../../src/problems.h"
#include "../../src/serialization.h"
#include "../../src/types.h"
#include "../exceptions.h"
#include "../utils.h"
#include "python_base.h"

using namespace boost::python;
using namespace pagmo;

// Wrapper to expose problems.
template <class Problem>
static inline class_<Problem,bases<problem::base> > problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def_pickle(generic_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	return retval;
}

BOOST_PYTHON_MODULE(_problem) {
	// Translate exceptions for this module.
	translate_exceptions();

	// Expose base problem class, including the virtual methods.
	typedef void (problem::base::*bounds_setter)(const decision_vector &);
	typedef void (problem::base::*bounds_setter_value)(const double &, const double &);
	typedef void (problem::base::*bounds_setter_vectors)(const decision_vector &, const decision_vector &);
	typedef constraint_vector (problem::base::*return_constraints)(const decision_vector &) const;
	typedef fitness_vector (problem::base::*return_fitness)(const decision_vector &) const;
	class_<problem::python_base>("_base",init<int,optional<int,int,int,int,const double &> >())
		.def(init<const decision_vector &, const decision_vector &, optional<int,int,int,int, const double &> >())
		.def(init<const problem::base &>())
		.def("__repr__", &problem::base::human_readable)
		.def("is_thread_safe",&problem::base::is_thread_safe)
		// Dimensions.
		.add_property("dimension", &problem::base::get_dimension, "Global dimension.")
		.add_property("f_dimension", &problem::base::get_f_dimension, "Fitness dimension.")
		.add_property("i_dimension", &problem::base::get_i_dimension, "Integer dimension.")
		.add_property("c_dimension", &problem::base::get_c_dimension, "Global constraints dimension.")
		.add_property("ic_dimension", &problem::base::get_ic_dimension, "Inequality constraints dimension.")
		// Constraints tolerance.
		.add_property("c_tol", &problem::base::get_c_tol, "Tolerance used in constraints analysis.")
		// Bounds.
		.add_property("lb",make_function(&problem::base::get_lb,return_value_policy<copy_const_reference>()), bounds_setter(&problem::base::set_lb), "Lower bounds.")
		.add_property("ub",make_function(&problem::base::get_ub,return_value_policy<copy_const_reference>()), bounds_setter(&problem::base::set_ub), "Upper bounds.")
		.def("set_bounds",bounds_setter_value(&problem::base::set_bounds),"Set all bounds to the input values.")
		.def("set_bounds",bounds_setter_vectors(&problem::base::set_bounds),"Set bounds to the input vectors.")
		.add_property("diameter",&problem::base::get_diameter, "Problem's diameter.")
		// Useful operators.
		.def(self == self)
		.def(self != self)
		.def("is_compatible",&problem::base::is_compatible,"Check compatibility with other problem.")
		// Comparisons.
		.def("compare_x",&problem::base::compare_x,"Compare decision vectors.")
		.def("verify_x",&problem::base::verify_x,"Check if decision vector is compatible with problem.")
		.def("compare_fc",&problem::base::compare_fc,"Simultaneous fitness-constraint comparison.")
		// Constraints.
		.def("compare_constraints",&problem::base::compare_constraints,"Compare constraint vectors.")
		.def("compute_constraints",return_constraints(&problem::base::compute_constraints),"Compute and return constraint vector.")
		.def("feasibility_x",&problem::base::feasibility_x,"Determine feasibility of decision vector.")
		.def("feasibility_c",&problem::base::feasibility_c,"Determine feasibility of constraint vector.")
		// Fitness.
		.def("objfun",return_fitness(&problem::base::objfun),"Compute and return fitness vector.")
		.def("compare_fitness",&problem::base::compare_fitness,"Compare fitness vectors.")
		// Virtual methods that can be (re)implemented.
		.def("__copy__", pure_virtual(&problem::base::clone))
		.def("get_name",&problem::base::get_name,&problem::python_base::default_get_name)
		.def("human_readable_extra", &problem::base::human_readable_extra, &problem::python_base::default_human_readable_extra)
		.def("_get_typename",&problem::python_base::get_typename)
		.def("_objfun_impl",&problem::python_base::py_objfun)
		.def("_equality_operator_extra",&problem::python_base::py_equality_operator_extra)
		.def("_compute_constraints_impl",&problem::python_base::py_compute_constraints_impl)
		.def_pickle(python_class_pickle_suite<problem::python_base>());
	
	// Ackley problem.
	problem_wrapper<problem::ackley>("ackley","Ackley function.")
		.def(init<int>());

	// Griewank problem.
	problem_wrapper<problem::griewank>("griewank","Griewank function.")
		.def(init<int>());
	
	// De Jong's problem.
	problem_wrapper<problem::dejong>("dejong","De Jong's function.")
		.def(init<int>());
	
	// michalewicz's problem.
	problem_wrapper<problem::michalewicz>("michalewicz","Michalewicz's function.")
		.def(init<int, optional<int> >());

	// GTOC1 problem.
	problem_wrapper<problem::gtoc_1>("gtoc_1","GTOC1 problem.");

	// GTOC2 problem.
	problem_wrapper<problem::gtoc_2>("gtoc_2","GTOC problem.")
		.def(init<int,int,int,int,optional<int,problem::gtoc_2::objective> >());

	// GTOC2's objectives enum.
	enum_<problem::gtoc_2::objective>("gtoc2_objective")
		.value("MASS",problem::gtoc_2::MASS)
		.value("TIME",problem::gtoc_2::TIME)
		.value("MASS_TIME",problem::gtoc_2::MASS_TIME);

	// Inventory problem.
	problem_wrapper<problem::inventory>("inventory","Inventory problem.")
		.def(init<int, int>());
	
	// Laplace problem.
	problem_wrapper<problem::laplace>("laplace","Laplace problem.")
		.def(init< const std::vector<int> &>());
	
	// Lennard Jones problem.
	problem_wrapper<problem::lennard_jones>("lennard_jones","Lennard Jones problem.")
		.def(init<int>());
	
	// Levy5 problem.
	problem_wrapper<problem::levy5>("levy5","Levy5 problem.")
		.def(init<int>());
	
	// Paraboloid problem.
	problem_wrapper<problem::paraboloid>("paraboloid","Multi-dimensional paraboloid miminisation.")
		.def(init<const decision_vector &, const decision_vector &>());

	// Rosenbrock problem.
	problem_wrapper<problem::rosenbrock>("rosenbrock","Multi-dimensional Rosenbrock function.")
		.def(init<int>());
	
	// Rosetta problem.
	problem_wrapper<problem::rosetta>("rosetta","Rosetta problem.");
	
	// Sagas problem.
	problem_wrapper<problem::sagas>("sagas","Sagas problem.");
	
	// Rastrigin problem.
	problem_wrapper<problem::rastrigin>("rastrigin","Generalised Rastrigin function.")
		.def(init<int>());

	// Golomb ruler problem.
	problem_wrapper<problem::golomb_ruler>("golomb_ruler","Optimal Golomb ruler search. Constructor from order of the Golomb ruler and maximum distance between consecutive marks.")
		.def(init<int,int>());

	// Schwefel function.
	problem_wrapper<problem::schwefel>("schwefel","Two-dimensional Schwefel function.")
		.def(init<int>());

	// Knapsack problem.
	problem_wrapper<problem::knapsack>("knapsack","Classical 01 knapsack problem. Constructor from vector of values, vector of weights and maximum weight.")
		.def(init<const std::vector<double> &, const std::vector<double> &, const double &>());

	// Himmelblau's function.
	problem_wrapper<problem::himmelblau>("himmelblau","Himmelblau's function.");

	// SNOPT toy problem.
	problem_wrapper<problem::snopt_toyprob>("snopt_toyprob","SNOPT toy problem.");

	// Branin's rcos funciotn.
	problem_wrapper<problem::branin>("branin","Branin's rcos function.");

	// Luksan Vlcek problem 1.
	problem_wrapper<problem::luksan_vlcek_1>("luksan_vlcek_1","Luksan Vlcek problem 1.")
		.def(init<int, optional<const double &, const double &> >());

	// Luksan Vlcek problem 2.
	problem_wrapper<problem::luksan_vlcek_2>("luksan_vlcek_2","Luksan Vlcek problem 2.")
		.def(init<int, optional<const double &, const double &> >());

	// Luksan Vlcek problem 3.
	problem_wrapper<problem::luksan_vlcek_3>("luksan_vlcek_3","Luksan Vlcek problem 3.")
		.def(init<int, optional<const double &, const double &> >());

	// String matching problem.
	problem_wrapper<problem::string_match>("string_match","String matching problem.")
		.def(init<const std::string &>());

	// Cassini 1.
	problem_wrapper<problem::cassini_1>("cassini_1","Cassini 1 interplanetary trajectory problem.");

	// Messenger full.
	problem_wrapper<problem::messenger_full>("messenger_full","Full Messenger problem.");

	// Cassini 2.
	problem_wrapper<problem::cassini_2>("cassini_2","Cassini 2 interplanetary trajectory problem.");
	
	// Traveling salesman problem
	problem_wrapper<problem::tsp>("tsp","Traveling salesman problem")
		.def(init<const std::vector<std::vector<double> > &>());
	
	// Tandem.
	problem_wrapper<problem::tandem>("tandem","Tandem problem.")
		.def(init< optional<int, double> >());
	
	// SCH
	problem_wrapper<problem::sch>("sch","Shaffer's study problem.")
		.def(init<>());
	// FON
	problem_wrapper<problem::fon>("fon","Fonseca and Fleming's study problem.")
		.def(init<>());
	// POL
	problem_wrapper<problem::pol>("pol","Poloni's study problem.")
		.def(init<>());
	// KUR
	problem_wrapper<problem::kur>("kur","Kursawe's study problem.")
		.def(init<>());
	// ZDT1
	problem_wrapper<problem::zdt1>("zdt1","ZDT1")
		.def(init<>());
	// ZDT2
	problem_wrapper<problem::zdt2>("zdt2","ZDT2")
		.def(init<>());
	// ZDT3
	problem_wrapper<problem::zdt3>("zdt1","ZDT3")
		.def(init<>());
	// ZDT4
	problem_wrapper<problem::zdt4>("zdt4","ZDT4")
		.def(init<>());
	// ZDT6
	problem_wrapper<problem::zdt6>("zdt6","ZDT6")
		.def(init<>());

	// Human mission to asteroids.
	problem_wrapper<problem::sample_return>("sample_return","Asteroid sample return problem.")
		.def(init<const ::kep_toolbox::planet &, optional<const double &> >())
		.def("get_delta_v",&problem::sample_return::get_delta_v);

	// Earth-planet problem.
	problem_wrapper<problem::earth_planet>("earth_planet","Earth-planet low-thrust problem.")
		.def(init< optional<int, std::string, const double &> >());

	// GTOC5 launch.
	problem_wrapper<problem::gtoc5_launch>("gtoc5_launch","GTOC5 launch phase.")
		.def(init< optional<int, int, problem::gtoc5_launch::objective, const double &> >());

	enum_<problem::gtoc5_launch::objective>("gtoc5_launch_objective")
		.value("MASS",problem::gtoc5_launch::MASS)
		.value("TIME",problem::gtoc5_launch::TIME);

	// GTOC5 randez-vouz.
	problem_wrapper<problem::gtoc5_rendezvous>("gtoc5_rendezvous","GTOC5 rendezvous phase.")
		.def(init< optional<int, int, int, const double &, const double &, const double &> >());

	enum_<problem::gtoc5_flyby::objective>("gtoc5_flyby_objective")
		.value("MASS",problem::gtoc5_flyby::MASS)
		.value("FINAL_EPOCH",problem::gtoc5_flyby::FINAL_EPOCH)
		.value("TIME",problem::gtoc5_flyby::TIME);

	// GTOC5 flyby.
	problem_wrapper<problem::gtoc5_flyby>("gtoc5_flyby","GTOC5 flyby phase.")
		.def(init< optional<int, int, int, int, const double &, const double &, problem::gtoc5_flyby::objective, const double &, const double &> >());

	// GTOC5 self flyby.
	problem_wrapper<problem::gtoc5_self_flyby>("gtoc5_self_flyby","GTOC5 self flyby phase.")
		.def(init< optional<int, int, const double &, const double &, const double &> >());

	// Register to_python conversion from smart pointer.
	register_ptr_to_python<problem::base_ptr>();
}
