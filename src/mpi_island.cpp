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

#include <string>

#include "algorithm/base.h"
#include "base_island.h"
#include "mpi_island.h"
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
mpi_island::mpi_island(const problem::base &p, const algorithm::base &a, int n, const double &migr_prob,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy, int processor):
	base_island(p,a,n,migr_prob,s_policy,r_policy)
{
	processor_id = processor;
}

/// Copy constructor.
/**
 * @see pagmo::base_island constructors.
 */
mpi_island::mpi_island(const mpi_island &isl):base_island(isl)
{
	processor_id = isl.processor_id;
}

/// Constructor from population.
/**
 * @see pagmo::base_island constructors.
 */
mpi_island::mpi_island(const population &pop, const algorithm::base &a, const double &migr_prob,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy, int processor):
	base_island(pop,a,migr_prob,s_policy,r_policy)
{
	processor_id = processor;
}

/// Assignment operator.
mpi_island &mpi_island::operator=(const mpi_island &isl)
{
	base_island::operator=(isl);
	return *this;
}

base_island_ptr mpi_island::clone() const
{
	return base_island_ptr(new mpi_island(*this));
}

// Thread safety attribute implementation.
bool mpi_island::is_thread_blocking() const
{
	return (m_pop.problem().is_blocking() || m_algo->is_blocking());
}

void mpi_island::perform_evolution(const algorithm::base &algo, population &pop) const
{
	algo.evolve(pop);
}

std::string mpi_island::get_name() const
{
	return "Local thread mpi_island";
}

}
