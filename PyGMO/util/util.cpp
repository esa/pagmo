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

#include"../../src/util/discrepancy.h"
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
}
