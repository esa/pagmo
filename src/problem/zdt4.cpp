/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#include <cmath>
#include <boost/math/constants/constants.hpp>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "base.h"
#include "zdt4.h"

namespace pagmo { namespace problem {

/**
 * Will construct ZDT4.
 * 
 * @param[in] dim integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
zdt4::zdt4(size_type dim):base_unc_mo(dim,0,2)
{
	// Set bounds.
	set_lb(-5.0);
	set_ub(5.0);
	set_lb(0,0.0);
	set_ub(0,1.0);
}

/// Clone method.
base_ptr zdt4::clone() const
{
	return base_ptr(new zdt4(*this));
}

/// Convergence metric for a decision_vector (0 = converged to the optimal front)
double zdt4::convergence_metric(const decision_vector &x) const
{
	double c = 0.0;
	double g = 0.0;

	for(problem::base::size_type j = 1; j < x.size(); ++j) {
		g += x[j]*x[j] - 10 * cos(4 * boost::math::constants::pi<double>() * x[j]);
	}
	c += 1 + 10 * (x.size()-1) + g;
	return c  - 1;
}

// ZDT4': lambda x: 1 + 10 * (len(x) - 1) + sum([xi**2 - 10*np.cos(4*np.pi*xi) for xi in x[1:]])

/// Implementation of the objective function.
void zdt4::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 2);
    pagmo_assert(x.size() == get_dimension());

	double g = 1 + 10 * (x.size() - 1);
	
	f[0] = x[0];

	for(problem::base::size_type i = 1; i < x.size(); ++i) {
		g += x[i]*x[i] - 10 * cos(4 * boost::math::constants::pi<double>() * x[i]);
	}
	
	f[1] = g * ( 1 - sqrt(x[0]/g) );
	
}

std::string zdt4::get_name() const
{
	return "ZDT4";
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::zdt4);
