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
#include "identity.h"

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n-dimensional and n-objective Identity problem.
 *
 * @param[in] n integer dimension of the problem with n objectives.
 *
 * @see problem::base constructors.
 */
identity::identity(int n):base(n, 0, n) {
	set_lb(-100);
	set_ub(100);
}

/// Clone method.
base_ptr identity::clone() const {
	return base_ptr(new identity(*this));
}

/// Implementation of the objective function.
void identity::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(x.size() == f.size());
	decision_vector::size_type n = x.size();

	for (decision_vector::size_type i=0; i<n; i++){
		f[i] = x[i];
	}
}

std::string identity::get_name() const
{
	return "Identity";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::identity)
