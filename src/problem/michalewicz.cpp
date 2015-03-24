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

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "michalewicz.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional Michalewicz problem.
 *
 * @param[in] n integer dimension of the problem.
 * @param[in] m sin factor exponent 
 *
 * @see problem::base constructors.
 */
michalewicz::michalewicz(int n, int m):base(n), m_m(m) {
	// Set bounds.
	set_lb(0);
	set_ub(boost::math::constants::pi<double>()); //pi
}

/// Clone method.
base_ptr michalewicz::clone() const {
	return base_ptr(new michalewicz(*this));
}

/// Implementation of the objective function.
void michalewicz::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	decision_vector::size_type n = x.size();
	double retval = 0.0;

	for (decision_vector::size_type i=0; i<n; i++){
		retval -= sin(x[i]) * pow(sin((i+1)*x[i]*x[i]/boost::math::constants::pi<double>()) , 2*m_m);
	}
	f[0] = retval;
}

std::string michalewicz::get_name() const
{
	return "Michalewicz";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::michalewicz)
