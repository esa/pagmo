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

#include <boost/mpi.hpp>
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
mpi_island::mpi_island(const problem::base &p, const algorithm::base &a, int n, int processor, const double &migr_prob,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
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
mpi_island::mpi_island(const population &pop, const algorithm::base &a, int processor, const double &migr_prob,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
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

// Method that perform the actual evolution for the island population, and is used to distribute the computation load over multiple processors
void mpi_island::perform_evolution(const algorithm::base &algo, population &pop) const
{
	boost::mpi::communicator world;
	// In case an island has been assigned an inexistant processor we throw an error
	if (processor_id >= world.size()) {
		pagmo_throw(value_error,"invalid processor assigned to islands");
	}
	// in case the processor is the root processor then we send the population that needs to evolve to the corresponding island
	if (world.rank() == 0) {
		world.send(processor_id,1,pop);
		std::cout << "sent island " << processor_id << " to process " << world.rank() << "." << std::endl;
	}
	// The processor assined to the island now receives the population from the root algorithm and runs the evolution and sends back the evolved island populations to the root processor
	if (world.rank() == processor_id) {
		std::cout << "evolving island " << processor_id << " on process " << world.rank() << "." << std::endl;
		world.recv(0,1,pop);
		algo.evolve(pop);		
		world.send(0,2,pop);
	}
	// the root processor now receives the evolved island populations and moves on
	if (world.rank() == 0) {
		std::cout << "received evolved population from island on process " << processor_id  << std::endl;
		world.recv(processor_id,2,pop);
	}

}

std::string mpi_island::get_name() const
{
	return "Local thread mpi_island";
}

}
