/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
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

#ifndef PAGMO_EXCEPTIONS_H
#define PAGMO_EXCEPTIONS_H

#include <exception>
#include <iostream>
#include <string>

#define _PAGMO_QUOTEME(x) #x
#define PAGMO_QUOTEME(x) _PAGMO_QUOTEME(x)
#define PAGMO_EXCTOR(s) (__FILE__ "," PAGMO_QUOTEME(__LINE__) ": " s ".")
#define pagmo_throw(ex,s) (throw ex(PAGMO_EXCTOR(s)))

class base_exception: public std::exception
{
	public:
		base_exception(const std::string &s): m_what(s) {
#ifndef PAGMO_BUILD_PYGMO
			std::cout << m_what << '\n';
#endif
		}
		virtual const char *what() const throw() {
			return m_what.c_str();
		}
		virtual ~base_exception() throw() {}
	protected:
		std::string m_what;
};

struct index_error: public base_exception {
	index_error(const std::string &s): base_exception(s) {}
};

struct value_error: public base_exception {
	value_error(const std::string &s): base_exception(s) {}
};

struct type_error: public base_exception {
	type_error(const std::string &s): base_exception(s) {}
};

struct assertion_error: public base_exception {
	assertion_error(const std::string &s): base_exception(s) {}
};

struct runtime_error: public base_exception {
	runtime_error(const std::string &s): base_exception(s) {}
};

#endif
