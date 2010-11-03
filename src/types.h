/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

// 30/01/10 Created by Francesco Biscani.

#ifndef PAGMO_TYPES_H
#define PAGMO_TYPES_H

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <iostream>
#include <vector>

#include "config.h"

namespace pagmo
{
	/// Decision vector type.
	typedef std::vector<double> decision_vector;
	/// Fitness vector type.
	typedef std::vector<double> fitness_vector;
	/// Constraint vector type.
	typedef std::vector<double> constraint_vector;

	/// Double to int converter.
	/**
	 * Will round a double to the nearest integer using "round half to even" for tie-breaking. Usage example:
	 @verbatim int i = double_to_int::convert(1.23); @endverbatim
	 * @see http://www.boost.org/doc/libs/release/libs/numeric/conversion/doc/html/index.html
	 * @see http://en.wikipedia.org/wiki/Rounding
	 */
	typedef boost::numeric::converter<int,double,boost::numeric::conversion_traits<int,double>,
		boost::numeric::def_overflow_handler,boost::numeric::RoundEven<double> > double_to_int;
}

// Give the possibility to disable stream overloads with appropriate #define.
#ifndef PAGMO_NO_STD_VECTOR_STREAM_OVERLOADS

namespace std
{
	/// Overload stream insertion operator for std::vector<double>.
	inline ostream &operator<<(ostream &os, const vector<double> &v)
	{
		os << '[';
		for (std::vector<double>::size_type i = 0; i < v.size(); ++i) {
			os << boost::lexical_cast<std::string>(v[i]);
			if (i != v.size() - 1) {
				os << ", ";
			}
		}
		os << ']';
		return os;
	}

	/// Overload stream insertion operator for std::vector<int>.
	inline ostream &operator<<(ostream &os, const vector<int> &v)
	{
		os << '[';
		for (std::vector<int>::size_type i = 0; i < v.size(); ++i) {
			os << v[i];
			if (i != v.size() - 1) {
				os << ", ";
			}
		}
		os << ']';
		return os;
	}

	/// Overload stream insertion operator for std::vector<unsigned int>.
	inline ostream &operator<<(ostream &os, const vector<unsigned int> &v)
	{
		os << '[';
		for (std::vector<unsigned int>::size_type i = 0; i < v.size(); ++i) {
			os << v[i];
			if (i != v.size() - 1) {
				os << ", ";
			}
		}
		os << ']';
		return os;
	}

	/// Overload stream insertion operator for std::vector<unsigned long int>.
	inline ostream &operator<<(ostream &os, const vector<unsigned long int> &v)
	{
		os << '[';
		for (std::vector<unsigned long int>::size_type i = 0; i < v.size(); ++i) {
			os << v[i];
			if (i != v.size() - 1) {
				os << ", ";
			}
		}
		os << ']';
		return os;
	}

	/// Overload stream insertion operator for std::vector<unsigned long long int>.
	inline ostream &operator<<(ostream &os, const vector<unsigned long long int> &v)
	{
		os << '[';
		for (std::vector<unsigned long long int>::size_type i = 0; i < v.size(); ++i) {
			os << v[i];
			if (i != v.size() - 1) {
				os << ", ";
			}
		}
		os << ']';
		return os;
	}
}

#endif

#endif
