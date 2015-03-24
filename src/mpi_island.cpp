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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <list>
#include <set>
#include <stdexcept>
#include <string>

#include "algorithm/base.h"
#include "base_island.h"
#include "exceptions.h"
#include "mpi_environment.h"
#include "mpi_island.h"
#include "migration/base_r_policy.h"
#include "migration/base_s_policy.h"
#include "population.h"
#include "problem/base.h"

namespace pagmo
{

boost::mutex mpi_island::m_proc_mutex;
boost::condition_variable mpi_island::m_proc_cond;
boost::mutex mpi_island::m_mpi_mutex;
std::list<mpi_island const *> mpi_island::m_queue;
boost::scoped_ptr<std::set<int> > mpi_island::m_available_processors;

/// Constructor from problem::base, algorithm::base, number of individuals, migration probability and selection/replacement policies.
/**
 * @see pagmo::base_island constructors.
 */
mpi_island::mpi_island(const algorithm::base &a, const problem::base &p, int n,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	base_island(a,p,n,s_policy,r_policy)
{}

/// Constructor from population.
/**
 * @see pagmo::base_island constructors.
 */
mpi_island::mpi_island(const algorithm::base &a, const population &pop,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	base_island(a,pop,s_policy,r_policy)
{}

/// Copy constructor.
/**
 * @see pagmo::base_island constructors.
 */
mpi_island::mpi_island(const mpi_island &isl):base_island(isl)
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

// Method that perform the actual evolution for the island population, and is used to distribute the computation load over multiple processors
void mpi_island::perform_evolution(const algorithm::base &algo, population &pop) const
{
	// Create copy of data to be transmitted - will use a std::pair for packing everything in a single object.
	const boost::shared_ptr<population> pop_copy(new population(pop));
	const algorithm::base_ptr algo_copy = algo.clone();
	const std::pair<boost::shared_ptr<population>,algorithm::base_ptr> out(pop_copy,algo_copy);
	const int processor = acquire_processor();
	if (mpi_environment::is_multithread()) {
		mpi_environment::send(out,processor);
	} else {
		// Lock down MPI mutex before sending.
		boost::lock_guard<boost::mutex> lock(m_mpi_mutex);
		mpi_environment::send(out,processor);
	}
	bool successful = false;
	boost::shared_ptr<population> in;
	try {
		if (mpi_environment::is_multithread()) {
			mpi_environment::recv(in,processor);
		} else {
			while (true) {
				{
					boost::lock_guard<boost::mutex> lock(m_mpi_mutex);
					if (mpi_environment::iprobe(processor)) {
						mpi_environment::recv(in,processor);
						break;
					}
				}
				boost::this_thread::sleep(boost::posix_time::milliseconds(10));
			}
		}
		successful = true;
	} catch (const boost::archive::archive_exception &e) {
		std::cout << "MPI Recv Error during island evolution using " << algo.get_name() << ": " << e.what() << std::endl;
	} catch (...) {
		std::cout << "MPI Recv Error during island evolution using " << algo.get_name() << ", unknown exception caught. :(" << std::endl;
	}
	release_processor(processor);
	if (successful) {
		// NOTE: implement via population::swap (to be written) in order to avoid
		// extra copying?
		pop = *in;
	}
}

/// Return a string identifying the island's type.
/**
 * @return the string "MPI island".
 */
std::string mpi_island::get_name() const
{
	return "MPI island";
}

void mpi_island::init_processors()
{
	// Make sure we are not called with the mutex not locked.
	pagmo_assert(!m_proc_mutex.try_lock());
	if (!m_available_processors) {
		m_available_processors.reset(new std::set<int>());
		pagmo_assert(mpi_environment::get_size() >= 2);
		// Fill in the available processors (the root processor is excluded).
		for (int i = 1; i < mpi_environment::get_size(); ++i) {
			m_available_processors->insert(i);
		}
	}
}

int mpi_island::acquire_processor() const
{
	// Lock down before doing anything else.
	boost::unique_lock<boost::mutex> lock(m_proc_mutex);
	init_processors();
	if (m_available_processors->empty() || !m_queue.empty()) {
		// Put self at the end of the queue.
		m_queue.push_back(this);
		while (*m_queue.begin() != this || m_available_processors->empty()) {
			m_proc_cond.wait(lock);
		}
		// this reached the head of the queue and there are available processors: pop it
		// from the head and proceed.
		m_queue.pop_front();
	}
	pagmo_assert(!m_available_processors->empty());
	const int retval = *m_available_processors->begin();
	m_available_processors->erase(m_available_processors->begin());
	return retval;
}

void mpi_island::release_processor(int n) const
{
	{
		boost::lock_guard<boost::mutex> lock(m_proc_mutex);
		init_processors();
		if (n <= 0 || n >= mpi_environment::get_size()) {
			pagmo_throw(std::runtime_error,"invalid processor id: the value is either non-positive or exceeding the size of the MPI world");
		}
		if (m_available_processors->find(n) != m_available_processors->end()) {
			pagmo_throw(std::runtime_error,"trying to release a processor which was never acquired");
		}
		// Re-insert the processor in the pool of available processors.
		m_available_processors->insert(n);
	}
	m_proc_cond.notify_all();
}

}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::mpi_island)
