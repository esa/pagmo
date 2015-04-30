/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

// Workaround for http://mail.python.org/pipermail/new-bugs-announce/2011-March/010395.html
#ifdef _WIN32
#include <cmath>
#endif
 
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
#include "../utils.h"
#include "python_base.h"
#include "python_base_stochastic.h"

#ifdef PAGMO_ENABLE_KEP_TOOLBOX
		#include <keplerian_toolbox/planet/base.h>
		#include "../../src/problem/tsp_ds.h"
#endif

using namespace boost::python;
using namespace pagmo;


// Transforms an Eigen Matrix into a std::vector<std::vector<double> >
// Used for the rotated meta-problem
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
// Used for the decomposed meta-problem
static inline fitness_vector compute_decomposed_fitness_wrapper(const problem::decompose& p, const fitness_vector &original_fit, const fitness_vector &weights) {
	fitness_vector retval(1);
	p.compute_decomposed_fitness(retval,original_fit,weights);
	return retval;
}

// Wrapper to expose stochastic meta-problems. (should we not here also derive from bases<problem::base_meta> ?)
template <class Problem>
static inline class_<Problem,bases<problem::base>,bases<problem::base_stochastic> > stochastic_problem_wrapper(const char *name, const char *descr)
{
	class_<Problem,bases<problem::base>,bases<problem::base_stochastic> > retval(name,descr,init<const Problem &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Problem>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Problem>);
	retval.def_pickle(python_class_pickle_suite<Problem>());
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
	retval.def_pickle(python_class_pickle_suite<Problem>());
	retval.def("cpp_loads", &py_cpp_loads<Problem>);
	retval.def("cpp_dumps", &py_cpp_dumps<Problem>);
	return retval;
}


BOOST_PYTHON_MODULE(_problem_meta) {
	common_module_init();

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

	// Decompose enums
	enum_<problem::decompose::method_type>("_decomposition_method")
		.value("WEIGHTED", problem::decompose::WEIGHTED)
		.value("TCHEBYCHEFF", problem::decompose::TCHEBYCHEFF)
		.value("BI", problem::decompose::BI);
	// Decomposition meta-problem
	meta_problem_wrapper<problem::decompose>("decompose","Decomposed problem")
		.def(init<const problem::base &, optional<problem::decompose::method_type, const std::vector<double> &, const std::vector<double> &, const bool> >())
		.def("compute_decomposed_fitness", &compute_decomposed_fitness_wrapper, 
		"Computes the fitness of the decomposed problem\n\n"
		"  USAGE:: w = prob.compute_decomposed_fitness(fit,weight)\n"
		"   - fit: multi-dimensional fitness\n"
		"   - weight: decomposition weights")
		.add_property("weights", make_function(&problem::decompose::get_weights, return_value_policy<copy_const_reference>()))
		.add_property("ideal_point", &problem::decompose::get_ideal_point, &problem::decompose::set_ideal_point,"the (z) point used to compute tchebycheff and bi decompositions");

	// Exposing enums of problem::noisy
	enum_<problem::noisy::noise_type>("_noise_distribution")
		.value("NORMAL", problem::noisy::NORMAL)
		.value("UNIFORM", problem::noisy::UNIFORM);
	// Noisy meta-problem
	stochastic_problem_wrapper<problem::noisy>("noisy", "Noisy problem")
		.def(init<const problem::base &, unsigned int, const double, const double, problem::noisy::noise_type, unsigned int>())
		.add_property("noise_param_first", &problem::noisy::get_param_first)
		.add_property("noise_param_second", &problem::noisy::get_param_second);

	// Robust meta-problem
	stochastic_problem_wrapper<problem::robust>("robust", "Robust problem")
		.def(init<const problem::base &,unsigned int, const double, unsigned int>())
		.add_property("rho", &problem::robust::get_rho);
}
