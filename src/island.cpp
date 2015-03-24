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

#include <string>

#include "algorithm/base.h"
#include "base_island.h"
#include "island.h"
#include "migration/base_r_policy.h"
#include "migration/base_s_policy.h"
#include "population.h"
#include "problem/base.h"

namespace pagmo
{

/// Constructor from problem::base, algorithm::base, number of individuals, migration probability and selection/replacement policies.
/**
 * @see pagmo::base_island constructors.
 */
island::island(const algorithm::base &a, const problem::base &p, int n,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	base_island(a,p,n,s_policy,r_policy)
{}

/// Copy constructor.
/**
 * @see pagmo::base_island constructors.
 */
island::island(const island &isl):base_island(isl)
{}

/// Constructor from population.
/**
 * @see pagmo::base_island constructors.
 */
island::island(const algorithm::base &a, const population &pop,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	base_island(a,pop,s_policy,r_policy)
{}

/// Assignment operator.
island &island::operator=(const island &isl)
{
	base_island::operator=(isl);
	return *this;
}

base_island_ptr island::clone() const
{
	return base_island_ptr(new island(*this));
}

void island::perform_evolution(const algorithm::base &algo, population &pop) const
{
	algo.evolve(pop);
}

std::string island::get_name() const
{
	return "Local thread island";
}

}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::island)
