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

// 13/02/2008: Initial version by Francesco Biscani.

#ifndef PYGMO_UTILS_H
#define PYGMO_UTILS_H

#include <Python.h>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/python/class.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/tuple.hpp>
#include <csignal>
#include <sstream>
#include <string>

#include "exceptions.h"

template <class T>
inline T Py_copy_from_ctor(const T &x)
{
	return T(x);
}

template <class T>
inline T Py_deepcopy_from_ctor(const T &x, boost::python::dict)
{
	return T(x);
}

// Serialization for python wrapper.
namespace boost { namespace serialization {

template <class Archive, class T>
void serialize(Archive &, boost::python::wrapper<T> &, const unsigned int)
{}

}}

// Generic pickle suite for C++ classes with default constructor extensible from Python.
// Difference from above is that we need to take care of handling the derived class' dict.
template <class T>
struct python_class_pickle_suite: boost::python::pickle_suite
{
	static boost::python::tuple getinitargs(const T &)
	{
		return boost::python::make_tuple();
	}
	static boost::python::tuple getstate(boost::python::object obj)
	{
		T const &x = boost::python::extract<T const &>(obj)();
		std::stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa << x;
		return boost::python::make_tuple(obj.attr("__dict__"),ss.str());
	}
	static void setstate(boost::python::object obj, boost::python::tuple state)
	{
		using namespace boost::python;
		T &x = extract<T &>(obj)();
		if (len(state) != 2)
		{
			PyErr_SetObject(PyExc_ValueError,("expected 2-item tuple in call to __setstate__; got %s" % state).ptr());
			throw_error_already_set();
		}
		// Restore the object's __dict__.
		boost::python::dict d = extract<boost::python::dict>(obj.attr("__dict__"))();
		d.update(state[0]);
		// Restore the internal state of the C++ object.
		const std::string str = extract<std::string>(state[1]);
		std::stringstream ss(str);
		boost::archive::text_iarchive ia(ss);
		ia >> x;
	}
	static bool getstate_manages_dict()
	{
		return true;
	}
};

template <class T>
inline void py_cpp_loads(T &x, const std::string &s)
{
	std::stringstream ss(s);
	boost::archive::text_iarchive ia(ss);
	ia >> x;
}

template <class T>
inline std::string py_cpp_dumps(const T &x)
{
	std::stringstream ss;
	boost::archive::text_oarchive oa(ss);
	oa << x;
	return ss.str();
}


#define common_module_init() \
/* Initialise Python thread support. */ \
PyEval_InitThreads(); \
/* Translate exceptions for this module. */ \
translate_exceptions(); \
/* Disable docstring C++ signature. */ \
boost::python::docstring_options local_docstring_options(true,true,false); 
#endif
