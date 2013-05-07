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

#ifndef P_EXCEPTIONS_H
#define P_EXCEPTIONS_H

#include <exception>
#include <iostream>
#include <string>

#define _P_EXCEPTION_QUOTEME(x) #x
#define P_EXCEPTION_QUOTEME(x) _P_EXCEPTION_QUOTEME(x)
#define P_EXCEPTION_EXCTOR(s) ((std::string(__FILE__ "," P_EXCEPTION_QUOTEME(__LINE__) ": ") + s) + ".")
#define P_EX_THROW(ex,s) (throw ex(P_EXCEPTION_EXCTOR(s)))

class p_base_exception: public std::exception {
	public:
		p_base_exception(const std::string &s):m_what(s) {}
		virtual const char *what() const throw() {
			return m_what.c_str();
		}
		virtual ~p_base_exception() throw() {}
	protected:
		std::string m_what;
};

struct index_error: public p_base_exception {
	index_error(const std::string &s): p_base_exception(s) {}
};

struct value_error: public p_base_exception {
	value_error(const std::string &s): p_base_exception(s) {}
};

struct io_error: public p_base_exception {
	io_error(const std::string &s): p_base_exception(s) {}
};

struct type_error: public p_base_exception {
	type_error(const std::string &s): p_base_exception(s) {}
};

struct assertion_error: public p_base_exception {
	assertion_error(const std::string &s): p_base_exception(s) {}
};

struct zero_division_error: public p_base_exception {
	zero_division_error(const std::string &s): p_base_exception(s) {}
};

struct not_implemented_error: public p_base_exception {
	not_implemented_error(const std::string &s): p_base_exception(s) {}
};

struct memory_error: public p_base_exception {
	memory_error(const std::string &s): p_base_exception(s) {}
};

#define P_EX_ASSERT(expr) \
if (!(expr)) { \
	P_EX_THROW(assertion_error,"assertion error"); \
}

#endif
