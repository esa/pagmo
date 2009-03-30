/***************************************************************************
 *   Copyright (C) 2009 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef BOOST_PYTHON_P_EXCEPTIONS_H
#define BOOST_PYTHON_P_EXCEPTIONS_H

#include <boost/python/exception_translator.hpp>

#include "p_exceptions.h"

inline void ie_translator(const index_error &ie)
{
	PyErr_SetString(PyExc_IndexError, ie.what());
}

inline void ve_translator(const value_error &ve)
{
	PyErr_SetString(PyExc_ValueError, ve.what());
}

inline void te_translator(const type_error &te)
{
	PyErr_SetString(PyExc_TypeError, te.what());
}

inline void ae_translator(const assertion_error &ae)
{
	PyErr_SetString(PyExc_AssertionError, ae.what());
}

inline void re_translator(const runtime_error &re)
{
	PyErr_SetString(PyExc_RuntimeError, re.what());
}

// Translate our C++ exceptions into Python exceptions.
inline void translate_p_exceptions()
{
	boost::python::register_exception_translator<index_error>(ie_translator);
	boost::python::register_exception_translator<value_error>(ve_translator);
	boost::python::register_exception_translator<type_error>(te_translator);
	boost::python::register_exception_translator<assertion_error>(ae_translator);
	boost::python::register_exception_translator<runtime_error>(re_translator);
}

#endif
