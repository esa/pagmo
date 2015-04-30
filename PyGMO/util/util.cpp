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
#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../../src/util/hypervolume.h"
#include "../../src/util/discrepancy.h"
#include "../../src/util/race_pop.h"
#include "../../src/util/race_algo.h"
#include "../utils.h"

#include "../../src/util/hv_algorithm/base.h"
#include "../../src/util/hv_algorithm/hv2d.h"
#include "../../src/util/hv_algorithm/hv3d.h"
#include "../../src/util/hv_algorithm/hv4d.h"
#include "../../src/util/hv_algorithm/wfg.h"
#include "../../src/util/hv_algorithm/bf_approx.h"
#include "../../src/util/hv_algorithm/bf_fpras.h"
#include "../../src/util/hv_algorithm/hoy.h"
#include "../../src/util/hv_algorithm/fpl.h"

using namespace boost::python;
using namespace pagmo;
using namespace pagmo::util;

namespace pagmo {namespace util{ namespace discrepancy{


class __PAGMO_VISIBLE py_simplex
{
	public:
		py_simplex(unsigned int dim, unsigned int count) : m_original_class(dim,count) {}
		std::vector<double> operator ()() {return m_original_class();}
		std::vector<double> operator ()(unsigned int n) {return m_original_class(n);}
	private:
		pagmo::util::discrepancy::simplex m_original_class;
};

class __PAGMO_VISIBLE py_sobol
{
	public:
		py_sobol(unsigned int dim, unsigned int count) : m_original_class(dim,count) {}
		std::vector<double> operator ()() {return m_original_class();}
		std::vector<double> operator ()(unsigned int n) {return m_original_class(n);}
	private:
		pagmo::util::discrepancy::sobol m_original_class;
};

class __PAGMO_VISIBLE py_lhs
{
	public:
		py_lhs(unsigned int dim, unsigned int count) : m_original_class(dim,count) {}
		std::vector<double> operator ()() {return m_original_class();}
		std::vector<double> operator ()(unsigned int n) {return m_original_class(n);}
	private:
		pagmo::util::discrepancy::lhs m_original_class;
};

class __PAGMO_VISIBLE py_halton
{
	public:
		py_halton(unsigned int dim, unsigned int count) : m_original_class(dim,count) {}
		std::vector<double> operator ()() {return m_original_class();}
		std::vector<double> operator ()(unsigned int n) {return m_original_class(n);}
	private:
		pagmo::util::discrepancy::halton m_original_class;
};

class __PAGMO_VISIBLE py_faure
{
	public:
		py_faure(unsigned int dim, unsigned int count) : m_original_class(dim,count) {}
		std::vector<double> operator ()() {return m_original_class();}
		std::vector<double> operator ()(unsigned int n) {return m_original_class(n);}
	private:
		pagmo::util::discrepancy::faure m_original_class;
};

}}}

template <class HVAlgorithm>
static inline class_<HVAlgorithm,bases<util::hv_algorithm::base> > algorithm_wrapper(const char *name, const char *descr)
{

	class_<HVAlgorithm,bases<util::hv_algorithm::base> > retval(name,descr,init<const HVAlgorithm &>());
	retval.def(init<>());
	return retval;
}

void expose_hv_algorithm()
{
	class_<util::hv_algorithm::base,boost::noncopyable>("_base",no_init)
		.def("get_name", &util::hv_algorithm::base::get_name);
	algorithm_wrapper<util::hv_algorithm::hv2d>("hv2d","hv2d algorithm.");
	algorithm_wrapper<util::hv_algorithm::hv3d>("hv3d","hv3d algorithm.");
	algorithm_wrapper<util::hv_algorithm::hv4d>("hv4d","hv4d algorithm.");
	algorithm_wrapper<util::hv_algorithm::fpl>("fpl","FPL algorithm.");
	algorithm_wrapper<util::hv_algorithm::hoy>("hoy","HOY algorithm.");
	class_<util::hv_algorithm::wfg, bases<util::hv_algorithm::base> >("wfg","WFG algorithm.", init<const unsigned int>());
	class_<util::hv_algorithm::bf_approx, bases<util::hv_algorithm::base> >("bf_approx","Bringmann-Friedrich approximated algorithm.", 
			init<const bool, const unsigned int, const double, const double, const double, const double, const double, const double>());
	class_<util::hv_algorithm::bf_fpras, bases<util::hv_algorithm::base> >("bf_fpras","Hypervolume approximation based on FPRAS", init<const double, const double>());
}

void expose_hypervolume()
{
	typedef double (util::hypervolume::*compute_custom)(const fitness_vector &, const util::hv_algorithm::base_ptr) const;
	typedef double (util::hypervolume::*compute_dynamic)(const fitness_vector &) const;

	typedef double (util::hypervolume::*exclusive_custom)(const unsigned int, const fitness_vector &, const util::hv_algorithm::base_ptr) const;
	typedef double (util::hypervolume::*exclusive_dynamic)(const unsigned int, const fitness_vector &) const;

	typedef unsigned int (util::hypervolume::*least_contributor_custom)(const fitness_vector &, const util::hv_algorithm::base_ptr) const;
	typedef unsigned int (util::hypervolume::*least_contributor_dynamic)(const fitness_vector &) const;

	typedef unsigned int (util::hypervolume::*greatest_contributor_custom)(const fitness_vector &, const util::hv_algorithm::base_ptr) const;
	typedef unsigned int (util::hypervolume::*greatest_contributor_dynamic)(const fitness_vector &) const;

	typedef std::vector<double> (util::hypervolume::*contributions_custom)(const fitness_vector &, const util::hv_algorithm::base_ptr) const;
	typedef std::vector<double> (util::hypervolume::*contributions_dynamic)(const fitness_vector &) const;

	class_<util::hypervolume>("hypervolume","Hypervolume class.", init<const std::vector<std::vector<double> > &, const bool >())
		.def(init<boost::shared_ptr<population>, const bool>())
		.def("compute", compute_custom(&util::hypervolume::compute), "Computes the hypervolume using the provided hypervolume algorithm.")
		.def("compute", compute_dynamic(&util::hypervolume::compute), "Computes the hypervolume.")
		.def("exclusive", exclusive_custom(&util::hypervolume::exclusive), "Computes the exclusive hypervolume using the provided hypervolume algorithm.")
		.def("exclusive", exclusive_dynamic(&util::hypervolume::exclusive), "Computes the exclusive hypervolume.")
		.def("least_contributor", least_contributor_custom(&util::hypervolume::least_contributor), "Get the least contributor of the hypervolume using provided hypervolume algorithm.")
		.def("least_contributor", least_contributor_dynamic(&util::hypervolume::least_contributor), "Get the least contributor of the hypervolume.")
		.def("greatest_contributor", greatest_contributor_custom(&util::hypervolume::greatest_contributor), "Get the greatest contributor of the hypervolume using provided hypervolume algorithm.")
		.def("greatest_contributor", greatest_contributor_dynamic(&util::hypervolume::greatest_contributor), "Get the greatest contributor of the hypervolume.")
		.def("contributions", contributions_custom(&util::hypervolume::contributions), "Get the contributions to the hypervolume by each point using provided hypervolume algorithm..")
		.def("contributions", contributions_dynamic(&util::hypervolume::contributions), "Get the contributions to the hypervolume by each point.")
		.def("get_nadir_point", &util::hypervolume::get_nadir_point)
		.def("set_copy_points", &util::hypervolume::set_copy_points)
		.def("get_copy_points", &util::hypervolume::get_copy_points)
		.def("get_points", &util::hypervolume::get_points)
		.def("set_verify", &util::hypervolume::set_verify)
		.def("get_verify", &util::hypervolume::get_verify);
}

// Main method containing all the juice of race_pop
static inline boost::python::tuple race_pop_run_return_tuple(
	racing::race_pop& race_obj,
	const pagmo::population::size_type n_final,
	const unsigned int min_trials,
	const unsigned int max_count,
	double delta,
	const std::vector<pagmo::population::size_type> &active_set,
	racing::race_pop::termination_condition term_cond,
	const bool race_best,
	const bool screen_output)
{
	std::pair<std::vector<pagmo::population::size_type>, unsigned int> res = race_obj.run(n_final, min_trials, max_count, delta, active_set, term_cond, race_best, screen_output);
	return boost::python::make_tuple(res.first, res.second);
}

// Main method containing all the juice of race_algo
static inline boost::python::tuple race_algo_run_return_tuple(
	racing::race_algo& race_obj,
	const pagmo::population::size_type n_final,
	const unsigned int min_trials,
	const unsigned int max_count,
	double delta,
	const std::vector<int> &active_set,
	const bool race_best,
	const bool screen_output)
{
	// boost python does not like vector of unsigned int......
	std::vector<unsigned int> active_set_(active_set.size());
	for(unsigned int i = 0; i < active_set.size(); i++){
		active_set_[i] = active_set[i];
	}
	std::pair<std::vector<unsigned int>, unsigned int> res = race_obj.run(n_final, min_trials, max_count, delta, active_set_, race_best, screen_output);
	// boost python does not like vector of unsigned int......
	std::vector<int> winners(res.first.size());
	for(unsigned int i = 0; i < res.first.size(); i++){
		winners[i] = res.first[i];
	}
	return boost::python::make_tuple(winners, res.second);
}

BOOST_PYTHON_MODULE(_util) {

	common_module_init();

	typedef std::vector<double> (discrepancy::py_lhs::*my_first_overload_l)() ;
	typedef std::vector<double> (discrepancy::py_lhs::*my_second_overload_l)(unsigned int) ;
	class_<discrepancy::py_lhs>("lhs", init<unsigned int , unsigned int>())
		.def("next", my_first_overload_l(&discrepancy::py_lhs::operator()))
		.def("next", my_second_overload_l(&discrepancy::py_lhs::operator()));

	typedef std::vector<double> (discrepancy::py_sobol::*my_first_overload_s)() ;
	typedef std::vector<double> (discrepancy::py_sobol::*my_second_overload_s)(unsigned int) ;
	class_<discrepancy::py_sobol>("sobol", init<unsigned int , unsigned int>())
		.def("next", my_first_overload_s(&discrepancy::py_sobol::operator()))
		.def("next", my_second_overload_s(&discrepancy::py_sobol::operator()));

	typedef std::vector<double> (discrepancy::py_simplex::*my_first_overload)() ;
	typedef std::vector<double> (discrepancy::py_simplex::*my_second_overload)(unsigned int) ;
	class_<discrepancy::py_simplex>("simplex", init<unsigned int , unsigned int>())
		.def("next", my_first_overload(&discrepancy::py_simplex::operator()))
		.def("next", my_second_overload(&discrepancy::py_simplex::operator()));
		
	typedef std::vector<double> (discrepancy::py_halton::*my_first_overload_h)() ;
	typedef std::vector<double> (discrepancy::py_halton::*my_second_overload_h)(unsigned int) ;
	class_<discrepancy::py_halton>("halton", init<unsigned int , unsigned int>())
		.def("next", my_first_overload_h(&discrepancy::py_halton::operator()))
		.def("next", my_second_overload_h(&discrepancy::py_halton::operator()));

	typedef std::vector<double> (discrepancy::py_faure::*my_first_overload_f)() ;
	typedef std::vector<double> (discrepancy::py_faure::*my_second_overload_f)(unsigned int) ;
	class_<discrepancy::py_faure>("faure", init<unsigned int , unsigned int>())
		.def("next", my_first_overload_f(&discrepancy::py_faure::operator()))
		.def("next", my_second_overload_f(&discrepancy::py_faure::operator()));

	// Racing
	enum_<racing::race_pop::termination_condition>("_termination_condition")
		.value("MAX_BUDGET", racing::race_pop::MAX_BUDGET)
		.value("MAX_DATA_COUNT", racing::race_pop::MAX_DATA_COUNT);
	class_<racing::race_pop>("race_pop", init<pagmo::population, unsigned int>())
		.def(init<unsigned int>())
		.def("run", &race_pop_run_return_tuple, "Race the individuals")
		.def("size", &racing::race_pop::size, "Returns number of individuals")
		.def("reset_cache", &racing::race_pop::reset_cache, "Clears the cache")
		.def("register_pop", &racing::race_pop::register_population, "Load a population into the race environment")
		.def("inherit_memory", &racing::race_pop::inherit_memory, "Transfer memory of identical decision vectors")
		.def("get_mean_fitness", &racing::race_pop::get_mean_fitness, "Returns the mean fitness of the individuals resulted from previously run race")
		.def("set_seed", &racing::race_pop::set_seed, "Set the ground seed of the race");

	// Required by race_algo
	//class_<std::vector<pagmo::algorithm::base_ptr> >("vector_of_algorithm_base_ptr")
	//	.def(vector_indexing_suite<std::vector<pagmo::algorithm::base_ptr>, true>());

	//class_<std::vector<pagmo::problem::base_ptr> >("vector_of_problem_base_ptr")
	//	.def(vector_indexing_suite<std::vector<pagmo::problem::base_ptr>, true>());

	class_<racing::race_algo>("race_algo", init<const std::vector<pagmo::algorithm::base_ptr> &, const pagmo::problem::base &, unsigned int, unsigned int>())
	.def(init<const std::vector<pagmo::algorithm::base_ptr> &, const std::vector<pagmo::problem::base_ptr> &, unsigned int, unsigned int>())
	.def("run", &race_algo_run_return_tuple, "Race the algorithms");
	
	// Hypervolumes
	expose_hypervolume();
	scope current;
	std::string submoduleName(extract<const char*>(current.attr("__name__")));
	submoduleName.append(".hv_algorithm");

	object submodule(borrowed(PyImport_AddModule(submoduleName.c_str())));
	current.attr("hv_algorithm") = submodule;
	scope submoduleScope = submodule;
	expose_hv_algorithm();
}
