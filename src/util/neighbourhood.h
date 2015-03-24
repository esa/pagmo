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

#ifndef PAGMO_UTIL_NEIGHBOURHOOD_H
#define PAGMO_UTIL_NEIGHBOURHOOD_H

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <boost/shared_ptr.hpp>

#include "../config.h"
#include "../population.h"
#include "../exceptions.h"

namespace pagmo{ namespace util {

/**
 * Utilities to build a neighbourhood graph given some vectors
 *
*/
namespace neighbourhood {

template<class T> class __PAGMO_VISIBLE sorter {
public:
	sorter(const std::vector<T> &v) : values(v) {}
	bool operator()(int a, int b) { return values[a] < values[b]; }
private:
	const std::vector<T> &values;
};

/// Sort according the the values in the values vector but return the permutation
template<class T> std::vector<population::size_type> order(const std::vector<T> &values)
{
	std::vector<pagmo::population::size_type> rv(values.size());
	int idx = 0;
	for (std::vector<pagmo::population::size_type>::iterator i = rv.begin(); i != rv.end(); i++)
		*i = idx++;
	std::sort(rv.begin(), rv.end(), sorter<T>(values));
	return rv;
}

/**
 * Build a neighbourhood graph for vectors of real numbers using the euclidian distance
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */
class __PAGMO_VISIBLE euclidian {
public:
	static void compute_neighbours(std::vector<std::vector<pagmo::population::size_type> > &, const std::vector<std::vector<double> > &);
	static double distance(const std::vector<double> &, const std::vector<double> &);
};

}}}

#endif
