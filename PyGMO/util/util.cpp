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
#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/utility.hpp>

#include "../../src/util/hypervolume.h"
#include"../../src/util/discrepancy.h"
#include "../utils.h"

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

void expose_hv_algorithm() {
	class_<util::hv_algorithm::base,boost::noncopyable>("_base",no_init)
		.def("get_name", &util::hv_algorithm::base::get_name);
	algorithm_wrapper<util::hv_algorithm::native2d>("native2d","Native2D algorithm.");
	algorithm_wrapper<util::hv_algorithm::beume3d>("beume3d","Beume3D algorithm.");
	algorithm_wrapper<util::hv_algorithm::hv4d>("hv4d","HV4D algorithm.");
	algorithm_wrapper<util::hv_algorithm::hoy>("hoy","HOY by Beume algorithm.");
	class_<util::hv_algorithm::wfg, bases<util::hv_algorithm::base> >("wfg","WFG algorithm.", init<const unsigned int>());
	class_<util::hv_algorithm::bf_approx, bases<util::hv_algorithm::base> >("bf_approx","Bringmann-Friedrich approximated algorithm.", 
			init<const bool, const unsigned int, const double, const double, const double, const double, const double, const double>());
	class_<util::hv_algorithm::bf_fpras, bases<util::hv_algorithm::base> >("bf_fpras","Hypervolume approximation based on FPRAS", init<const double, const double>());
}

void expose_hypervolume() {

	typedef double (util::hypervolume::*compute_custom)(const fitness_vector &, const util::hv_algorithm::base_ptr) const;
	typedef double (util::hypervolume::*compute_dynamic)(const fitness_vector &) const;

	typedef double (util::hypervolume::*exclusive_custom)(const unsigned int, const fitness_vector &, const util::hv_algorithm::base_ptr) const;
	typedef double (util::hypervolume::*exclusive_dynamic)(const unsigned int, const fitness_vector &) const;

	typedef unsigned int (util::hypervolume::*least_contributor_custom)(const fitness_vector &, const util::hv_algorithm::base_ptr) const;
	typedef unsigned int (util::hypervolume::*least_contributor_dynamic)(const fitness_vector &) const;

	typedef unsigned int (util::hypervolume::*greatest_contributor_custom)(const fitness_vector &, const util::hv_algorithm::base_ptr) const;
	typedef unsigned int (util::hypervolume::*greatest_contributor_dynamic)(const fitness_vector &) const;


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
		.def("get_nadir_point", &util::hypervolume::get_nadir_point)
		.def("set_copy_points", &util::hypervolume::set_copy_points)
		.def("get_copy_points", &util::hypervolume::get_copy_points)
		.def("set_verify", &util::hypervolume::set_verify)
		.def("get_verify", &util::hypervolume::get_verify);
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

	expose_hypervolume();

	scope current;
	std::string submoduleName(extract<const char*>(current.attr("__name__")));
	submoduleName.append(".hv_algorithm");

	object submodule(borrowed(PyImport_AddModule(submoduleName.c_str())));
	current.attr("hv_algorithm") = submodule;
	scope submoduleScope = submodule;
	expose_hv_algorithm();
}
