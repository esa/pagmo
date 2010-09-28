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

// 01/02/10 Created by Francesco Biscani.

#include <string>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "paraboloid.h"

namespace pagmo { namespace problem {

static const double paraboloid_default_lb[1] = {-1};
static const double paraboloid_default_ub[1] = {1};

/// Default constructor.
/**
 * Will construct a one-dimensional problem with bounds [-1,1].
 */
paraboloid::paraboloid():base(paraboloid_default_lb,paraboloid_default_ub) {}

/// Constructor from lower/upper bounds.
/**
 * @see problem::base
 */
paraboloid::paraboloid(const decision_vector &lb, const decision_vector &ub):base(lb,ub) {}

/// Clone method.
base_ptr paraboloid::clone() const
{
	return base_ptr(new paraboloid(*this));
}

/// Implementation of the objective function.
void paraboloid::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	typedef decision_vector::size_type size_type;
	const size_type size = x.size();
	f[0] = 0;
	for (size_type i = 0; i < size; ++i) {
		f[0] += x[i] * x[i];
	}
}

std::string paraboloid::get_name() const
{
	return "Paraboloid";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::paraboloid);
