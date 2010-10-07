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

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "algorithm/base.h"
#include "base_island.h"
#include "exceptions.h"
#include "mpi_island.h"
#include "migration/base_r_policy.h"
#include "migration/base_s_policy.h"
#include "population.h"
#include "problem/base.h"

namespace pagmo
{

boost::mutex mpi_island::m_mutex;
boost::scoped_ptr<std::set<int> > mpi_island::m_available_processors;

/// Constructor from problem::base, algorithm::base, number of individuals, migration probability and selection/replacement policies.
/**
 * @see pagmo::base_island constructors.
 */
mpi_island::mpi_island(const problem::base &p, const algorithm::base &a, int n, const double &migr_prob,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	base_island(p,a,n,migr_prob,s_policy,r_policy)
{}

/// Copy constructor.
/**
 * @see pagmo::base_island constructors.
 */
mpi_island::mpi_island(const mpi_island &isl):base_island(isl)
{}

/// Constructor from population.
/**
 * @see pagmo::base_island constructors.
 */
mpi_island::mpi_island(const population &pop, const algorithm::base &a, const double &migr_prob,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	base_island(pop,a,migr_prob,s_policy,r_policy)
{}

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

bool mpi_island::is_blocking_impl() const
{
	return false;
}

// Method that perform the actual evolution for the island population, and is used to distribute the computation load over multiple processors
void mpi_island::perform_evolution(const algorithm::base &algo, population &pop) const
{
	// Create copy of data to be transmitted - will use a std::pair for packing everything in a single object.
	const boost::shared_ptr<population> pop_copy(new population(pop));
	const algorithm::base_ptr algo_copy = algo.clone();
	const std::pair<const boost::shared_ptr<population>, const algorithm::base_ptr> out(pop_copy,algo_copy);
	MPI_Status status;
	int processor, flag, size;
	{
		std::stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa << out;
		// Copy the content of the stream to a char vector.
		const std::string buffer_str(ss.str());
		std::vector<char> buffer_char(buffer_str.begin(),buffer_str.end());
		// Cast needed to safely convert to an MPI data type.
		size = boost::numeric_cast<int>(buffer_char.size());
		// Lock down.
		lock_type lock(m_mutex);
		processor = acquire_processor();
// std::cout << "master sending size " << processor << '\n';
		MPI_Send(static_cast<void *>(&size),1,MPI_INT,processor,0,MPI_COMM_WORLD);
// std::cout << "master sent size " << processor << '\n';
// std::cout << "master sending payload " << processor << '\n';
		MPI_Send(static_cast<void *>(&buffer_char[0]),size,MPI_CHAR,processor,1,MPI_COMM_WORLD);
// std::cout << "master sent payload " << processor << '\n';
	}
	std::vector<char> buffer_char;
	while (true) {
		{
			lock_type lock(m_mutex);
// std::cout << "master quering " << processor << '\n';
			MPI_Iprobe(processor,0,MPI_COMM_WORLD,&flag,&status);
			if (flag) {
// std::cout << "master receiving size " << processor << '\n';
				MPI_Recv(static_cast<void *>(&size),1,MPI_INT,processor,0,MPI_COMM_WORLD,&status);
// std::cout << "master received size " << processor << '\n';
				// Prepare buffer.
				buffer_char.resize(boost::numeric_cast<std::vector<char>::size_type>(size),0);
// std::cout << "master receiving payload " << processor << '\n';
				MPI_Recv(static_cast<void *>(&buffer_char[0]),size,MPI_CHAR,processor,1,MPI_COMM_WORLD,&status);
// std::cout << "master received payload " << processor << '\n';
				release_processor(processor);
				break;
			}
		}
		boost::this_thread::sleep(boost::posix_time::milliseconds(100));
	}
	const std::string buffer_str(buffer_char.begin(),buffer_char.end());
	std::stringstream ss(buffer_str);
	boost::archive::text_iarchive ia(ss);
	std::pair<boost::shared_ptr<population>,algorithm::base_ptr> in;
	ia >> in;
	pop = *in.first;
}

std::string mpi_island::get_name() const
{
	return "MPI island";
}

void mpi_island::init_processors()
{
	if (!m_available_processors) {
		m_available_processors.reset(new std::set<int>());
		boost::mpi::communicator world;
		if (world.size() < 2) {
			pagmo_throw(std::runtime_error,"the size of the MPI world must be at least 2");
		}
		// Fill in the available processors (the root processor is excluded).
		for (int i = 1; i < world.size(); ++i) {
			m_available_processors->insert(i);
		}
	}
}

int mpi_island::acquire_processor()
{
	init_processors();
	if (m_available_processors->empty()) {
		pagmo_throw(std::runtime_error,"no more processors are available");
	}
	const int retval = *m_available_processors->begin();
	m_available_processors->erase(m_available_processors->begin());
	return retval;
}

void mpi_island::release_processor(int n)
{
	init_processors();
	boost::mpi::communicator world;
	if (n <= 0 || n >= world.size()) {
		pagmo_throw(std::runtime_error,"invalid processor id: the value is either non-positive or exceeding the size of the MPI world");
	}
	if (m_available_processors->find(n) != m_available_processors->end()) {
		pagmo_throw(std::runtime_error,"trying to release a processor which was never acquired");
	}
	// Re-insert the processor in the pool of available processors.
	m_available_processors->insert(n);
}

}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::mpi_island);
