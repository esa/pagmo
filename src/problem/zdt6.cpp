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
#include "base_unc_mo.h"
#include "zdt6.h"

namespace pagmo { namespace problem {

/**
 * Will construct ZDT6.
 *
 * @param[in] dim integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
zdt6::zdt6(size_type dim):base_unc_mo(dim,0,2)
{
	// Set bounds.
	set_lb(0.0);
	set_ub(1.0);
}

/// Clone method.
base_ptr zdt6::clone() const
{
	return base_ptr(new zdt6(*this));
}

/// Convergence metric for a decision_vector (0 = converged to the optimal front)
double zdt6::convergence_metric(const decision_vector &x) const
{
	double c = 0.0;
	double g = 0.0;

	for(problem::base::size_type j = 1; j < x.size(); ++j) {
		g += x[j];
	}
	c += 1 + 9 * pow((g / (x.size()-1)),0.25);

	return c - 1;
}


/// Implementation of the objective function.
void zdt6::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 2);
	pagmo_assert(x.size() == get_dimension());

	double g = 0;

	f[0] = 1 - exp(-4*x[0])*pow(sin(6*boost::math::constants::pi<double>()*x[0]),6);

	for(problem::base::size_type i = 1; i < x.size(); ++i) {
		g += x[i];
	}
	g = 1 + 9 * pow((g/(x.size()-1)),0.25);

	f[1] = g * ( 1 - (f[0]/g)*(f[0]/g));

}

std::string zdt6::get_name() const
{
	return "ZDT6";
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::zdt6);
