/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#include <Python.h>
#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>
#include <cstddef>
#include <string>

#include "../../src/exceptions.h"
#include "../../src/problems.h"
#include "../../src/serialization.h"
#include "../../src/types.h"
#include "../../src/config.h"
#include "../../src/problem/tsp_ds.h"
#include "../utils.h"
#include "python_base.h"
#include "python_base_stochastic.h"

#ifdef PAGMO_ENABLE_KEP_TOOLBOX
		#include "../../src/keplerian_toolbox/keplerian_toolbox.h"
#endif

using namespace boost::python;
using namespace pagmo;


// wrapper for the find_subsequence method of tsp_cs
static inline tuple find_subsequence_wrapper_cs(const problem::tsp_cs& p, const decision_vector& tour)
{
	double retval_p, retval_l;
	decision_vector::size_type retval_it_l, retval_it_r;
    p.find_subsequence(tour,retval_p,retval_l,retval_it_l,retval_it_r);
	return boost::python::make_tuple(retval_p,retval_l,retval_it_l,retval_it_r);
}

// wrapper for the find_subsequence method of tsp_ads
static inline tuple find_subsequence_wrapper_ds(const problem::tsp_ds& p, const decision_vector& tour, const bool static_computations)
{
	double retval_p, retval_l;
	decision_vector::size_type retval_it_l, retval_it_r;
    p.find_subsequence(tour,retval_p,retval_l,retval_it_l,retval_it_r, static_computations);
	return boost::python::make_tuple(retval_p,retval_l,retval_it_l,retval_it_r);
}

// Transforms an Eigen Matrix into a std::vector<std::vector<double> >
std::vector<std::vector<double> > get_rotation_matrix_from_eigen(const problem::rotated & p) {
	Eigen::MatrixXd rot = p.get_rotation_matrix();
	pagmo_assert(rot.cols()==rot.rows());
	size_t dim = rot.cols();
	std::vector<double> dummy(dim,0);
	std::vector<std::vector<double> > retval(dim,dummy);

	for (size_t i=0; i<dim;++i){
		for (size_t j=0; j<dim;++j){
			retval[i][j] =rot(i,j);
		}
	}
	return retval;
}

// wrapper of a decompose method
static inline fitness_vector compute_decomposed_fitness_wrapper(const problem::decompose& p, const fitness_vector &original_fit, const fitness_vector &weights) {
	fitness_vector retval(1);
	p.compute_decomposed_fitness(retval,original_fit,weights);
	return retval;
}

// Wrapper to expose problems.
template <class Problem>
static inline class_<Problem,bases<problem::base> > problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Problem>);
	retval.def_pickle(generic_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	return retval;
}

// Wrapper to expose stochastic problems.
template <class Problem>
static inline class_<Problem,bases<problem::base>,bases<problem::base_stochastic> > stochastic_problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base>,bases<problem::base_stochastic> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Problem>);
	retval.def_pickle(generic_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	retval.add_property("seed",&problem::base_stochastic::get_seed,&problem::base_stochastic::set_seed,
		"Random seed used in the objective function evaluation.");
	return retval;
}

// Wrapper to expose meta problems.
template <class Problem>
static inline class_<Problem,bases<problem::base>,bases<problem::base_meta> > meta_problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base>,bases<problem::base_meta> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Problem>);
	retval.def_pickle(generic_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	return retval;
}

// Wrapper to expose unconstrained multi-objective problems.
template <class Problem>
static inline class_<Problem,bases<problem::base>,bases<problem::base_unc_mo> > unc_mo_problem_wrapper(const char *name, const char *descr)
{
	// allows for overload of p_distance
	double (problem::base_unc_mo::*p_dist_o1)(const decision_vector &) const = &problem::base_unc_mo::p_distance;
	double (problem::base_unc_mo::*p_dist_o2)(const population &) const = &problem::base_unc_mo::p_distance;
	
	// actual class exposition
	class_<Problem,bases<problem::base>,bases<problem::base_unc_mo> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Problem>);
	retval.def_pickle(generic_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	retval.def("p_distance", p_dist_o1,
		"The p distance is a convergence metric measuring the distance of a population or individual from the pareto front.\n"
		"It is typically 0.0 if the individuals lie on the Pareto-front.\n\n" 
		"  USAGE: x = prob.p_distance(pop)\n"
		"  USAGE: x = prob.p_distance(x)\n\n"
		"  * pop: population to evaluate\n"
		"  * x: chromosome to evaluate");
	retval.def("p_distance", p_dist_o2);
	return retval;
}

// Wrapper to expose TSP problems.
template <class Problem>
static inline class_<Problem,bases<problem::base>,bases<problem::base_tsp> > tsp_problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base>,bases<problem::base_tsp> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Problem>);
	retval.def_pickle(generic_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	retval.def("full2cities", &problem::base_tsp::full2cities);
	retval.def("cities2full", &problem::base_tsp::cities2full);
	retval.def("randomkeys2cities", &problem::base_tsp::randomkeys2cities);
	retval.def("cities2randomkeys", &problem::base_tsp::cities2randomkeys);
	retval.add_property("encoding", &problem::base_tsp::get_encoding);
	retval.add_property("n_cities", &problem::base_tsp::get_n_cities);
	return retval;
}


BOOST_PYTHON_MODULE(_problem) {
	common_module_init();

	// Expose base problem class, including the virtual methods.
	typedef void (problem::base::*bounds_setter)(const decision_vector &);
	typedef void (problem::base::*bounds_setter_value)(const double &, const double &);
	typedef void (problem::base::*bounds_setter_vectors)(const decision_vector &, const decision_vector &);
	typedef void (problem::base::*best_x_setter)(const std::vector<decision_vector>&);
	typedef constraint_vector (problem::base::*return_constraints)(const decision_vector &) const;
	typedef fitness_vector (problem::base::*return_fitness)(const decision_vector &) const;
    class_<problem::python_base, boost::noncopyable>("_base",init<int,optional<int,int,int,int,const std::vector<double> &> >())
		.def(init<const decision_vector &, const decision_vector &, optional<int,int,int,int, const double &> >())
		.def(init<int,int,int,int,int,const double>())
		.def("__repr__", &problem::base::human_readable)
		// Dimensions.
		.add_property("dimension", &problem::base::get_dimension, "Global dimension.")
		.add_property("f_dimension", &problem::base::get_f_dimension, "Fitness dimension.")
		.add_property("i_dimension", &problem::base::get_i_dimension, "Integer dimension.")
		.add_property("c_dimension", &problem::base::get_c_dimension, "Global constraints dimension.")
		.add_property("ic_dimension", &problem::base::get_ic_dimension, "Inequality constraints dimension.")
		// Constraints tolerance.
		.add_property("c_tol", make_function(&problem::base::get_c_tol,return_value_policy<copy_const_reference>()), "Tolerance used in constraints analysis.")
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
		.def("reset_caches",&problem::base::reset_caches,"Resets the internal caching system of PyGMO that stores previos calls to the objective function/ constraint function. This method should be called whenever a problem object is changed and the change affects the objective function.")
		// Comparisons.
		.def("compare_x",&problem::base::compare_x,"Compare decision vectors.")
		.def("verify_x",&problem::base::verify_x,"Check if decision vector is compatible with problem.")
		.def("compare_fc",&problem::base::compare_fc,"Simultaneous fitness-constraint comparison.")
		// Constraints.
		.def("compare_constraints",&problem::base::compare_constraints,"Compare constraint vectors.")
		.def("compute_constraints",return_constraints(&problem::base::compute_constraints),"Compute and return constraint vector.")
		.def("test_constraint",&problem::base::test_constraint,"Determine feasibility of the i-th constraint.")
		.def("feasibility_x",&problem::base::feasibility_x,"Determine feasibility of decision vector.")
		.def("feasibility_c",&problem::base::feasibility_c,"Determine feasibility of constraint vector.")
		// Fitness.
		.def("objfun",return_fitness(&problem::base::objfun),"Compute and return fitness vector.")
		.def("compare_fitness",&problem::base::compare_fitness,"Compare fitness vectors.")
		// Virtual methods that can be (re)implemented.
		.def("get_name",&problem::base::get_name,&problem::python_base::default_get_name)
		.def("human_readable_extra", &problem::base::human_readable_extra, &problem::python_base::default_human_readable_extra)
		.def("_get_typename",&problem::python_base::get_typename)
		.def("_objfun_impl",&problem::python_base::py_objfun)
		.def("_equality_operator_extra",&problem::python_base::py_equality_operator_extra)
		.def("_compute_constraints_impl",&problem::python_base::py_compute_constraints_impl)
		.def("_compare_constraints_impl",&problem::python_base::py_compare_constraints_impl)
		.def("_compare_fc_impl",&problem::python_base::py_compare_fc_impl)
		.def("_compare_fitness_impl",&problem::python_base::py_compare_fitness_impl)
		// Best known solution
		.add_property("best_x",make_function(&problem::base::get_best_x,return_value_policy<copy_const_reference>()), best_x_setter(&problem::base::set_best_x), "Best known decision vector(s).")
		.add_property("best_f",make_function(&problem::base::get_best_f,return_value_policy<copy_const_reference>()),"Best known fitness vector(s).")
		.add_property("best_c",make_function(&problem::base::get_best_c,return_value_policy<copy_const_reference>()),"Best known constraints vector(s).")
		.add_property("fevals",&problem::base::get_fevals,"Number of function evaluations.")
		.add_property("cevals",&problem::base::get_cevals,"Number of constraints evaluations.")
		.def_pickle(python_class_pickle_suite<problem::python_base>());

	// Expose base stochastic problem class, including the virtual methods. Here we explicitly
	// tell python that these objects can be passed where problem::base is expected .... in the previous
	// case (problem::python_base) this was not necessary as problem::python_base derives from
	// boost::python::wrapper<base> which informs python on the correct inheritance.
	class_<problem::python_base_stochastic, boost::noncopyable, bases<problem::base> >("_base_stochastic",init<int,optional<unsigned int> >())
		.def("__repr__", &problem::base::human_readable)
		// Dimensions.
		.add_property("dimension", &problem::base::get_dimension, "Global dimension.")
		.add_property("f_dimension", &problem::base::get_f_dimension, "Fitness dimension.")
		.add_property("i_dimension", &problem::base::get_i_dimension, "Integer dimension.")
		.add_property("c_dimension", &problem::base::get_c_dimension, "Global constraints dimension.")
		.add_property("ic_dimension", &problem::base::get_ic_dimension, "Inequality constraints dimension.")
		// Constraints tolerance.
		.add_property("c_tol", make_function(&problem::base::get_c_tol,return_value_policy<copy_const_reference>()), "Tolerance used in constraints analysis.")
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
		.def("test_constraint",&problem::base::test_constraint,"Determine feasibility of the i-th constraint.")
		.def("feasibility_x",&problem::base::feasibility_x,"Determine feasibility of decision vector.")
		.def("feasibility_c",&problem::base::feasibility_c,"Determine feasibility of constraint vector.")
		// Fitness.
		.def("objfun",return_fitness(&problem::base::objfun),"Compute and return fitness vector.")
		.def("compare_fitness",&problem::base::compare_fitness,"Compare fitness vectors.")
		// Seed.
		.add_property("seed",&problem::base_stochastic::get_seed,&problem::base_stochastic::set_seed,"Random seed used in the objective function evaluation.")
		// Virtual methods that can be (re)implemented.
		.def("get_name",&problem::base::get_name,&problem::python_base_stochastic::default_get_name)
		.def("human_readable_extra", &problem::base::human_readable_extra, &problem::python_base_stochastic::default_human_readable_extra)
		.def("_get_typename",&problem::python_base_stochastic::get_typename)
		.def("_objfun_impl",&problem::python_base_stochastic::py_objfun)
		.def("_equality_operator_extra",&problem::python_base_stochastic::py_equality_operator_extra)
		.def("_compute_constraints_impl",&problem::python_base_stochastic::py_compute_constraints_impl)
		.def("_compare_constraints_impl",&problem::python_base_stochastic::py_compare_constraints_impl)
		.def("_compare_fc_impl",&problem::python_base_stochastic::py_compare_fc_impl)
		.def("_compare_fitness_impl",&problem::python_base_stochastic::py_compare_fitness_impl)
		// Best known solution
		.add_property("best_x",make_function(&problem::base::get_best_x,return_value_policy<copy_const_reference>()), best_x_setter(&problem::base::set_best_x), "Best known decision vector(s).")
		.add_property("best_f",make_function(&problem::base::get_best_f,return_value_policy<copy_const_reference>()),"Best known fitness vector(s).")
		.add_property("best_c",make_function(&problem::base::get_best_c,return_value_policy<copy_const_reference>()),"Best known constraints vector(s).")
		.def_pickle(python_class_pickle_suite<problem::python_base_stochastic>());

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

	// Identity problem.
	stochastic_problem_wrapper<problem::identity>("identity","Identity problem.")
		.def(init<int>());

	// Inventory problem.
	stochastic_problem_wrapper<problem::inventory>("inventory","Inventory problem.")
		.def(init<optional<int, int, int> >());

	// Lennard Jones problem.
	problem_wrapper<problem::lennard_jones>("lennard_jones","Lennard Jones problem.")
		.def(init<int>());

	// Levy5 problem.
	problem_wrapper<problem::levy5>("levy5","Levy5 problem.")
		.def(init<int>());

	// Rosenbrock problem.
	problem_wrapper<problem::rosenbrock>("rosenbrock","Multi-dimensional Rosenbrock function.")
		.def(init<int>());

	// Rastrigin problem.
	problem_wrapper<problem::rastrigin>("rastrigin","Generalised Rastrigin function.")
		.def(init<int>());

	// Golomb ruler problem.
	problem_wrapper<problem::golomb_ruler>("golomb_ruler","Optimal Golomb ruler search. Constructor from order of the Golomb ruler and maximum distance between consecutive marks.")
		.def(init<int,int>());

	// Schwefel function.
	problem_wrapper<problem::schwefel>("schwefel","Two-dimensional Schwefel function.")
		.def(init<int>());

	// Himmelblau's function.
	problem_wrapper<problem::himmelblau>("himmelblau","Himmelblau's function.");
	
	// Bukin's f6 unction.
	problem_wrapper<problem::bukin>("bukin","Bukin's f6 function.");

	// CEC2006 Competition Problems.
	problem_wrapper<problem::cec2006>("cec2006","CEC2006 Competition Problems.")
		.def(init<int>());

	// Pressure vessel design.
	problem_wrapper<problem::pressure_vessel>("pressure_vessel","Pressure vessel design problem.");

	// Weight minimization of a tension compression string.
	problem_wrapper<problem::tens_comp_string>("tens_comp_string","Weight minimization of a tension compression string problem.");

	// Welded beam design Problems.
	problem_wrapper<problem::welded_beam>("welded_beam","Welded beam design problem.");

	// CEC2009 Competition Problems.
	problem_wrapper<problem::cec2009>("cec2009","CEC2009 Competition Problems.")
		.def(init<int, problem::base::size_type, bool>());

	// CEC2013 Competition Problems.
	problem_wrapper<problem::cec2013>("cec2013","CEC2013 Competition Problems.")
		.def(init<unsigned int, problem::base::size_type, const std::string&>())
		.add_property("origin_shift", &problem::cec2013::origin_shift, "Returns the origin shift used to define the problem");

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
		.def(init<const std::string &>())
		.def("pretty", &problem::string_match::pretty);

	// Travelling salesman problem (TSP) encoding enums
	enum_<problem::base_tsp::encoding_type>("_tsp_encoding")
		.value("FULL", problem::base_tsp::FULL)
		.value("RANDOMKEYS", problem::base_tsp::RANDOMKEYS)
		.value("CITIES", problem::base_tsp::CITIES);

	// Travelling salesman problem (TSP)
	tsp_problem_wrapper<problem::tsp>("tsp","Travelling salesman problem (TSP and ATSP)")
		.def(init<const std::vector<std::vector<double> > &, const problem::base_tsp::encoding_type &>())
		.add_property("weights", make_function(&problem::tsp::get_weights, return_value_policy<copy_const_reference>()));

	// Travelling salesman problem, vehicle routing problem with limited capacity variant (TSP-VRPLC)
	tsp_problem_wrapper<problem::tsp_vrplc>("tsp_vrplc","Vehicle routing problem with limited capacity (TSP-VRPLC)")
		.def(init<const std::vector<std::vector<double> > &, const problem::base_tsp::encoding_type &, const double&>())
		.def("return_tours",&problem::tsp_vrplc::return_tours,"Compute and return list of tours.")
		.add_property("weights", make_function(&problem::tsp_vrplc::get_weights, return_value_policy<copy_const_reference>()))
		.add_property("capacity", make_function(&problem::tsp_vrplc::get_capacity, return_value_policy<copy_const_reference>()));

	// Travelling salesman problem, city-selection variant (TSP-CS)
	tsp_problem_wrapper<problem::tsp_cs>("tsp_cs","City-selection Travelling Salesman Problem (TSP-CS)")
		.def(init<const std::vector<std::vector<double> > &, const std::vector<double>&, const double, const problem::base_tsp::encoding_type &>())
		.def("find_city_subsequence", &find_subsequence_wrapper_cs)
		.add_property("weights", make_function(&problem::tsp_cs::get_weights, return_value_policy<copy_const_reference>()))
		.add_property("values",  make_function(&problem::tsp_cs::get_values, return_value_policy<copy_const_reference>()))
		.add_property("max_path_length",  &problem::tsp_cs::get_max_path_length);

	// SCH
	problem_wrapper<problem::sch>("sch","Shaffer's study problem.");
	// FON
	problem_wrapper<problem::fon>("fon","Fonseca and Fleming's study problem.");
	// POL
	problem_wrapper<problem::pol>("pol","Poloni's study problem.");
	// KUR
	problem_wrapper<problem::kur>("kur","Kursawe's study problem.")
		.def(init<int>());
	// DTLZ
	unc_mo_problem_wrapper<problem::dtlz>("dtlz","DTLZ benchmark problem suite.")
		.def(init<optional<size_t, const size_t, fitness_vector::size_type, const size_t> >());
	// ZDT
	unc_mo_problem_wrapper<problem::zdt>("zdt","ZDT")
		.def(init<optional<size_t, size_t> >());
	
	// Meta-problems

	
// Death penalty enums
	enum_<problem::death_penalty::method_type>("_death_method_type")
		.value("SIMPLE", problem::death_penalty::SIMPLE)
		.value("KURI", problem::death_penalty::KURI)
		.value("WEIGHTED", problem::death_penalty::WEIGHTED);

	// Death penalty meta-problem
	meta_problem_wrapper<problem::death_penalty>("death_penalty","Constrained death penalty problem")
		.def(init<optional<const problem::base &, problem::death_penalty::method_type, const std::vector<double> &> >());

	// con2mo penalty enums
	enum_<problem::con2mo::method_type>("_con2mo_method_type")
		.value("OBJ_CSTRS", problem::con2mo::OBJ_CSTRS)
		.value("OBJ_CSTRSVIO", problem::con2mo::OBJ_CSTRSVIO)
		.value("OBJ_EQVIO_INEQVIO", problem::con2mo::OBJ_EQVIO_INEQVIO);

	// Constrained to multi-objective meta-problem.
	meta_problem_wrapper<problem::con2mo>("con2mo","Constrained to multi-objective problem")
		.def(init<optional<const problem::base &, problem::con2mo::method_type> >());

	// con2uncon penalty enums
	enum_<problem::con2uncon::method_type>("_con2uncon_method_type")
		.value("OPTIMALITY", problem::con2uncon::OPTIMALITY)
		.value("FEASIBILITY", problem::con2uncon::FEASIBILITY);

	// Constrained to unconstrained meta-problem.
	meta_problem_wrapper<problem::con2uncon>("con2uncon","Constrained to unconstrained problem")
		.def(init<optional<const problem::base &, problem::con2uncon::method_type> >());

	// Shifted meta-problem
	meta_problem_wrapper<problem::shifted>("shifted","Shifted problem")
		.def(init<const problem::base &>())
		.def(init<const problem::base &, std::vector<double> >())
		.def(init<const problem::base &, double>())
		.add_property("shift_vector",make_function(&problem::shifted::get_shift_vector,return_value_policy<copy_const_reference>()))
		.def("deshift",&problem::shifted::deshift);
		
	// Scaled meta-problem
	meta_problem_wrapper<problem::scaled>("scaled","Scaled problem")
		.def(init<const problem::base &, fitness_vector >())
		.add_property("units",make_function(&problem::scaled::get_units,return_value_policy<copy_const_reference>()))
		.def("descale",&problem::scaled::descale);
		
	// Rotated meta-problem
	meta_problem_wrapper<problem::rotated>("rotated","Rotated problem")
		.def(init<const problem::base &>())
		.def(init<const problem::base &, Eigen::MatrixXd >())
		.add_property("rotation_matrix",&get_rotation_matrix_from_eigen)
		.def("derotate",&problem::rotated::derotate);
		
	// Normalized meta-problem
	meta_problem_wrapper<problem::normalized>("normalized","Normalized problem")
		.def(init<const problem::base &>())
		.def("denormalize", &problem::normalized::denormalize);

	// Decomposition meta-problem
	// Exposing enums of problem::decompose
	enum_<problem::decompose::method_type>("_decomposition_method")
		.value("WEIGHTED", problem::decompose::WEIGHTED)
		.value("TCHEBYCHEFF", problem::decompose::TCHEBYCHEFF)
		.value("BI", problem::decompose::BI);
	meta_problem_wrapper<problem::decompose>("decompose","Decomposed problem")
		.def(init<const problem::base &, optional<problem::decompose::method_type, const std::vector<double> &, const std::vector<double> &, const bool> >())
		.def("compute_decomposed_fitness", &compute_decomposed_fitness_wrapper, 
		"Computes the fitness of the decomposed problem\n\n"
		"  USAGE:: w = prob.compute_decomposed_fitness(fit,weight)\n"
		"   - fit: multi-dimensional fitness\n"
		"   - weight: decomposition weights")
		.add_property("weights", make_function(&problem::decompose::get_weights, return_value_policy<copy_const_reference>()))
		.add_property("ideal_point", &problem::decompose::get_ideal_point, &problem::decompose::set_ideal_point,"the (z) point used to compute tchebycheff and bi decompositions");


	// Noisy meta-problem
	// Exposing enums of problem::noisy
	enum_<problem::noisy::noise_type>("_noise_distribution")
		.value("NORMAL", problem::noisy::NORMAL)
		.value("UNIFORM", problem::noisy::UNIFORM);

	stochastic_problem_wrapper<problem::noisy>("noisy", "Noisy problem")
		.def(init<const problem::base &, unsigned int, const double, const double, problem::noisy::noise_type, unsigned int>())
		.add_property("noise_param_first", &problem::noisy::get_param_first)
		.add_property("noise_param_second", &problem::noisy::get_param_second);

	// Robust meta-problem
	stochastic_problem_wrapper<problem::robust>("robust", "Robust problem")
		.def(init<const problem::base &,unsigned int, const double, unsigned int>())
		.add_property("rho", &problem::robust::get_rho);


#ifdef PAGMO_ENABLE_KEP_TOOLBOX
	// Asteroid Sample Return (also used fot human missions to asteroids)
//	problem_wrapper<problem::sample_return>("sample_return","Asteroid sample return problem.")
//		.def(init<const ::kep_toolbox::planet &, optional<const double &> >())
//		.def("get_delta_v",&problem::sample_return::get_delta_v);

	// Earth-planet problem.
	//problem_wrapper<problem::earth_planet>("earth_planet","Earth-planet low-thrust problem.")
		//.def(init< optional<int, std::string, const double &> >());

	// GTOC5 launch.
	//problem_wrapper<problem::gtoc5_launch>("gtoc5_launch","GTOC5 launch phase.")
		//.def(init< optional<int, int, problem::gtoc5_launch::objective, const double &> >());

	//enum_<problem::gtoc5_launch::objective>("gtoc5_launch_objective")
		//.value("MASS",problem::gtoc5_launch::MASS)
		//.value("TIME",problem::gtoc5_launch::TIME);

	// GTOC5 randez-vouz.
	//problem_wrapper<problem::gtoc5_rendezvous>("gtoc5_rendezvous","GTOC5 rendezvous phase.")
		//.def(init< optional<int, int, int, const double &, const double &, const double &> >());

	//enum_<problem::gtoc5_flyby::objective>("gtoc5_flyby_objective")
		//.value("MASS",problem::gtoc5_flyby::MASS)
		//.value("FINAL_EPOCH",problem::gtoc5_flyby::FINAL_EPOCH)
		//.value("TIME",problem::gtoc5_flyby::TIME);

	// GTOC5 flyby.
	//problem_wrapper<problem::gtoc5_flyby>("gtoc5_flyby","GTOC5 flyby phase.")
		//.def(init< optional<int, int, int, int, const double &, const double &, problem::gtoc5_flyby::objective, const double &, const double &> >());

	// GTOC5 self flyby.
	//problem_wrapper<problem::gtoc5_self_flyby>("gtoc5_self_flyby","GTOC5 self flyby phase.")
		//.def(init< optional<int, int, const double &, const double &, const double &> >());

	// GTOC1 problem.
	problem_wrapper<problem::gtoc_1>("gtoc_1","GTOC 1 problem (chemical approximation).");

	// GTOC2 problem.
	problem_wrapper<problem::gtoc_2>("gtoc_2","GTOC 2 problem (LT model).")
		.def(init<optional<int,int,int,int,int,problem::gtoc_2::objective> >());

	// GTOC2's objectives enum.
	enum_<problem::gtoc_2::objective>("_gtoc_2_objective")
		.value("MASS",problem::gtoc_2::MASS)
		.value("TIME",problem::gtoc_2::TIME)
		.value("MASS_TIME",problem::gtoc_2::MASS_TIME);

	// Laplace problem.
	problem_wrapper<problem::laplace>("laplace","Laplace problem.")
		.def(init< const std::vector<int> &>());

	// Cassini 1.
	problem_wrapper<problem::cassini_1>("cassini_1","Cassini 1 interplanetary trajectory problem.")
		.def(init<optional<unsigned int> >());
	// Messenger full.
	problem_wrapper<problem::messenger_full>("messenger_full","Full Messenger problem.");

	// Cassini 2.
	problem_wrapper<problem::cassini_2>("cassini_2","Cassini 2 interplanetary trajectory problem.");

	// Rosetta problem.
	problem_wrapper<problem::rosetta>("rosetta","Rosetta problem.");

	// Sagas problem.
	problem_wrapper<problem::sagas>("sagas","Sagas problem.");

	// Tandem.
	problem_wrapper<problem::tandem>("tandem","Tandem problem.")
		.def(init< optional<int, double> >())
		.def("pretty", &problem::tandem::pretty);

	problem_wrapper<problem::mga_1dsm_alpha>("mga_1dsm_alpha", "A Multiple Gravity Assist with 1 Deep Space Manouvre problem")
		.def(init< optional<std::vector<kep_toolbox::planet_ptr>, kep_toolbox::epoch, kep_toolbox::epoch, double, double, double, double, bool, bool, bool> >())
		.def("pretty", &problem::mga_1dsm_alpha::pretty)
		.def("set_tof", &problem::mga_1dsm_alpha::set_tof)
		.def("get_tof", &problem::mga_1dsm_alpha::get_tof)
		.def("set_launch_window", &problem::mga_1dsm_alpha::set_launch_window)
		.def("set_vinf", &problem::mga_1dsm_alpha::set_vinf)
		.def("get_sequence", &problem::mga_1dsm_alpha::get_sequence);

	problem_wrapper<problem::mga_1dsm_tof>("mga_1dsm_tof", "A Multiple Gravity Assist with 1 Deep Space Manouvre problem")
		.def(init< optional<std::vector<kep_toolbox::planet_ptr>, kep_toolbox::epoch, kep_toolbox::epoch, std::vector<boost::array<double,2> >, double, double, bool, bool, bool> >())
		.def("pretty", &problem::mga_1dsm_tof::pretty, (arg("x"),arg("extended_output") = false))
		.def("set_tof", &problem::mga_1dsm_tof::set_tof)
		.def("get_tof", &problem::mga_1dsm_tof::get_tof)
		.def("set_launch_window", &problem::mga_1dsm_tof::set_launch_window)
		.def("set_vinf", &problem::mga_1dsm_tof::set_vinf)
		.def("get_sequence", &problem::mga_1dsm_tof::get_sequence);
		
	problem_wrapper<problem::mga_incipit>("mga_incipit", "Jupiter capture problem from the first part of gtoc6")
		.def(init< optional< std::vector<kep_toolbox::planet_ptr>, kep_toolbox::epoch, kep_toolbox::epoch, std::vector<std::vector<double> > > >())
		.def("pretty", &problem::mga_incipit::pretty)
		.def("get_sequence", &problem::mga_incipit::get_sequence)
		.add_property("tof",make_function(&problem::mga_incipit::get_tof, return_value_policy<copy_const_reference>()), &problem::mga_incipit::set_tof,"bound on the times of flight for the different legs");

	problem_wrapper<problem::mga_incipit_cstrs>("mga_incipit_cstrs", "Jupiter capture problem from the first part of gtoc6 (constrained version)")
		.def(init< optional< std::vector<kep_toolbox::planet_ptr>, kep_toolbox::epoch, kep_toolbox::epoch, std::vector<std::vector<double> >, const double, const double > >())
		.def("pretty", &problem::mga_incipit_cstrs::pretty)
		.def("get_sequence", &problem::mga_incipit_cstrs::get_sequence)
		.add_property("tof",make_function(&problem::mga_incipit_cstrs::get_tof, return_value_policy<copy_const_reference>()), &problem::mga_incipit_cstrs::set_tof,"bound on the times of flight for the different legs");

	problem_wrapper<problem::mga_part>("mga_part", "A part of the Jupiter moon tour from gtoc6")
		.def(init< optional <std::vector<kep_toolbox::planet_ptr>, std::vector<std::vector<double> >, kep_toolbox::epoch, kep_toolbox::array3D > >())
		.def("pretty", &problem::mga_part::pretty)
		.def("get_sequence", &problem::mga_part::get_sequence)
		.add_property("vinf_in",make_function(&problem::mga_part::get_vinf_in, return_value_policy<copy_const_reference>()), &problem::mga_part::set_vinf_in,"initial incoming relative spacecraft velocity")
		.add_property("t0",make_function(&problem::mga_part::get_t0, return_value_policy<copy_const_reference>()), &problem::mga_part::set_t0, "start epoch")
		.add_property("tof",make_function(&problem::mga_part::get_tof, return_value_policy<copy_const_reference>()), &problem::mga_part::set_tof,"bounds on the times of flight for the different legs")
		.add_property("betas",&problem::mga_part::get_betas, &problem::mga_part::set_betas,"bounds on the beta angles for the different legs")
		.add_property("rps",&problem::mga_part::get_rps, &problem::mga_part::set_rps,"bounds on the periplanet heights for the different legs");

	// Travelling salesman problem, Asteroids / Debris Selection TSP (TSP-ADS)
	tsp_problem_wrapper<problem::tsp_ds>("tsp_ds","Debris Selection TSP (TSP-DS)")
		.def(init<optional<const std::vector<kep_toolbox::planet_ptr>&, const std::vector<double>&, const double, const std::vector<double>&, const problem::base_tsp::encoding_type & > >())
		.def("find_subsequence", &find_subsequence_wrapper_ds)
        .add_property("planets", make_function(&problem::tsp_ds::get_planets, return_value_policy<copy_const_reference>()) )
        .add_property("values", make_function(&problem::tsp_ds::get_values, return_value_policy<copy_const_reference>()) )
        .add_property("epochs", make_function(&problem::tsp_ds::get_epochs, return_value_policy<copy_const_reference>()), "epoch schedule")
        .add_property("max_DV", &problem::tsp_ds::get_max_DV );
#endif

#ifdef PAGMO_ENABLE_GSL
	// Spheres Problems
	stochastic_problem_wrapper<problem::spheres>("mit_spheres", "Spheres problem, a neurocontroller for the MIT test-bed (absolute perception-action)")
		.def(init< optional<int,int,double,unsigned int, bool, double, std::vector<double> > >())
		.def("post_evaluate", &problem::spheres::post_evaluate)
		.def("simulate", &problem::spheres::simulate)
		.def("get_nn_weights", &problem::spheres::get_nn_weights)
		.add_property("seed",&problem::spheres::get_seed,&problem::spheres::set_seed,"Random seed used in the objective function evaluation.");

	//problem_wrapper<problem::spheres_q>("spheres_q", "Spheres problem, a neurocontroller for the MIT test-bed (body-axis perception-action)")
	//	.def(init< optional<int,int,double,unsigned int> >())
	//	.def("post_evaluate", &problem::spheres_q::post_evaluate)
	//	.def("simulate", &problem::spheres_q::simulate);
#endif

	// Register to_python conversion from smart pointer.
	register_ptr_to_python<problem::base_ptr>();
}
