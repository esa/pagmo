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
#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>
#include <cstddef>
#include <string>

#include "../../src/exceptions.h"
#include "../../src/keplerian_toolbox/keplerian_toolbox.h"
#include "../../src/problems.h"
#include "../../src/types.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace pagmo;

template <class Problem>
static inline class_<Problem,bases<problem::base> > problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base> > retval(name,descr,init<const Problem &>());
	return retval;
}

struct python_problem: problem::base, wrapper<problem::base>
{
	python_problem(int n, int ni = 0, int nf = 1, int nc = 0, int nic = 0, const double &c_tol = 0):
		problem::base(n,ni,nf,nc,nic,c_tol) {}
	python_problem(const decision_vector &lb, const decision_vector &ub, int ni = 0, int nf = 1, int nc = 0, int nic = 0, const double &c_tol = 0):
		problem::base(lb,ub,ni,nf,nc,nic,c_tol) {}
	python_problem(const problem::base &p):problem::base(p) {}
	problem::base_ptr clone() const
	{
		if (override f = this->get_override("__copy__")) {
			return f();
		}
		pagmo_throw(not_implemented_error,"problem cloning has not been implemented");
		return problem::base_ptr();
	}
	std::string get_name() const
	{
		if (override f = this->get_override("get_name")) {
			return f();
		}
		return problem::base::get_name();
	}
	std::string default_get_name() const
	{
		return this->problem::base::get_name();
	}
	void objfun_impl(fitness_vector &f, const decision_vector &x) const
	{
		f = py_objfun(x);
	}
	fitness_vector py_objfun(const decision_vector &x) const
	{
		if (override f = this->get_override("_objfun_impl")) {
			return f(x);
		}
		pagmo_throw(not_implemented_error,"objective function has not been implemented");
		return fitness_vector();
	}
	bool is_blocking() const
	{
		return true;
	}
	std::string get_typename() const
	{
		return this->get_override("_get_typename")();
	}
	bool py_equality_operator_extra(const problem::base &p) const
	{
		if (override f = this->get_override("_equality_operator_extra")) {
			return f(p);
		}
		return problem::base::equality_operator_extra(p);
	}
	bool equality_operator_extra(const base &p) const
	{
		// NOTE: here the dynamic cast is safe because in base equality we already checked the C++ type.
		if (get_typename() != dynamic_cast<const python_problem &>(p).get_typename()) {
			return false;
		}
		return py_equality_operator_extra(p);
	}
	std::string human_readable_extra() const
	{
		return py_human_readable_extra();
	}
	std::string py_human_readable_extra() const
	{
		if (override f = this->get_override("_human_readable_extra")) {
			return f();
		}
		return problem::base::human_readable_extra();
	}
	void compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
	{
		if (this->get_override("_compute_constraints_impl")) {
			// If the function is overridden, use it.
			c = py_compute_constraints_impl(x);
		} else {
			// Otherwise, avoid memory allocation by calling the base function directly.
			problem::base::compute_constraints_impl(c,x);
		}
	}
	constraint_vector py_compute_constraints_impl(const decision_vector &x) const
	{
		override f = this->get_override("_compute_constraints_impl");
		pagmo_assert(f);
		return f(x);
	}
};

BOOST_PYTHON_MODULE(_problem) {
	// Translate exceptions for this module.
	translate_exceptions();

	// Expose base problem class, including the virtual methods.
	typedef void (problem::base::*bounds_setter)(const decision_vector &);
	typedef void (problem::base::*bounds_setter_value)(const double &, const double &);
	typedef void (problem::base::*bounds_setter_vectors)(const decision_vector &, const decision_vector &);
	typedef constraint_vector (problem::base::*return_constraints)(const decision_vector &) const;
	typedef fitness_vector (problem::base::*return_fitness)(const decision_vector &) const;
	class_<python_problem>("_base",no_init)
		.def(init<int,optional<int,int,int,int,const double &> >())
		.def(init<const decision_vector &, const decision_vector &, optional<int,int,int,int, const double &> >())
		.def(init<const problem::base &>())
		.def("__repr__", &problem::base::human_readable)
		.def("__copy__", &problem::base::clone)
		.def("is_blocking",&problem::base::is_blocking)
		.def("_get_typename",&python_problem::get_typename)
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
		.def("get_name",&problem::base::get_name,&python_problem::default_get_name)
		.def("_objfun_impl",&python_problem::py_objfun)
		.def("_human_readable_extra",&python_problem::py_human_readable_extra)
		.def("_equality_operator_extra",&python_problem::py_equality_operator_extra)
		.def("_compute_constraints_impl",&python_problem::py_compute_constraints_impl);
	
	// Ackley problem.
	problem_wrapper<problem::ackley>("ackley","Ackley function.")
		.def(init<int>());

	// Griewank problem.
	problem_wrapper<problem::griewank>("griewank","Griewank function.")
		.def(init<int>());

	// GTOC1 problem.
	problem_wrapper<problem::gtoc_1>("gtoc_1","GTOC1 problem.")
		.def(init<>());

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
		.def(init<>())
		.def(init<const decision_vector &, const decision_vector &>());

	// Rosenbrock problem.
	problem_wrapper<problem::rosenbrock>("rosenbrock","Multi-dimensional Rosenbrock function.")
		.def(init<int>());
	
	// NSGA-II fon problem.
	problem_wrapper<problem::nsga_ii_fon>("nsga_ii_fon","NSGA-II FON problem.")
		.def(init<>());
	
	// NSGA-II sch problem.
	problem_wrapper<problem::nsga_ii_sch>("nsga_ii_sch","NSGA-II SCH problem.")
		.def(init<>());
	
	// Rosetta problem.
	problem_wrapper<problem::rosetta>("rosetta","Rosetta problem.")
		.def(init<>());
	
	// Sagas problem.
	problem_wrapper<problem::sagas>("sagas","Sagas problem.")
		.def(init<>());
	
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
	problem_wrapper<problem::himmelblau>("himmelblau","Himmelblau's function.")
		.def(init<>());

	// SNOPT toy problem.
	problem_wrapper<problem::snopt_toyprob>("snopt_toyprob","SNOPT toy problem.")
		.def(init<>());

	// Branin's rcos funciotn.
	problem_wrapper<problem::branin>("branin","Branin's rcos function.")
		.def(init<>());

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

	// String matching problem (multi-objective version).
	problem_wrapper<problem::string_match_mo>("string_match_mo","String matching problem (multi-objective version).")
		.def(init<const std::string &>());

	// Cassini 1.
	problem_wrapper<problem::cassini_1>("cassini_1","Cassini 1 interplanetary trajectory problem.")
		.def(init<>());

	// Messenger full.
	problem_wrapper<problem::messenger_full>("messenger_full","Full Messenger problem.")
		.def(init<>());

	// Cassini 2.
	problem_wrapper<problem::cassini_2>("cassini_2","Cassini 2 interplanetary trajectory problem.")
		.def(init<>());
	
	// Tandem.
	problem_wrapper<problem::tandem>("tandem","Tandem problem.")
		.def(init< optional<int, double> >());

	// Human mission to asteroids.
	problem_wrapper<problem::sample_return>("sample_return","Asteroid sample return problem.")
		.def(init<const ::kep_toolbox::planet &, optional<const double &> >())
		.def("get_delta_v",&problem::sample_return::get_delta_v);

	// Function for the total number of objective function evaluations.
	def("objfun_calls",&problem::objfun_calls,"Return the total number of calls to the objective function.");
	def("reset_objfun_calls",&problem::reset_objfun_calls,"Reset the total number of calls to the objective function.");

	// Register to_python conversion from smart pointer.
	register_ptr_to_python<problem::base_ptr>();
}
