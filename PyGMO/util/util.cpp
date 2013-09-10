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
#include<boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/tuple.hpp>

#include"../../src/util/discrepancy.h"
#include "../../src/util/race_pop.h"
#include "../utils.h"

using namespace boost::python;
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


// Main method containing all the juice
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

BOOST_PYTHON_MODULE(_util) {
	common_module_init();
	
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
}
