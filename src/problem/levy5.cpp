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

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <string>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "levy5.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional Levy problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
levy5::levy5(int n):base(n)
{
	if (n < 2) {
		pagmo_throw(value_error,"the Levy5 problem's dimension must be at least 2");
	}
	// Set bounds.
	set_lb(-100);
	set_ub(100);
}

/// Clone method.
base_ptr levy5::clone() const
{
	return base_ptr(new levy5(*this));
}

/// Implementation of the objective function.
void levy5::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	decision_vector::size_type n = x.size();
	double isum = 0.0;
	double jsum = 0.0;
	f[0] = 0;

	for ( decision_vector::size_type j=0; j<n; j+=2 ) {
		for ( int i=1; i<=5; i++ ) {
			isum += (double)(i) * cos((double)(i-1)*x[j] + (double)(i));
			jsum += (double)(i) * cos((double)(i+1)*x[j+1] + (double)(i));
		}
	}

	f[0] = isum*jsum;
	for ( decision_vector::size_type j=0; j<n; j+=2 )
		f[0] += pow(x[j] + 1.42513,2) + pow(x[j+1] + 0.80032,2);

}

std::string levy5::get_name() const
{
	return "Levy5";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::levy5)
