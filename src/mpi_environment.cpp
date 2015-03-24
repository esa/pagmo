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

#include <boost/shared_ptr.hpp>
#include <cstdlib>
#include <mpi.h>
#include <stdexcept>
#include <utility>

#include "exceptions.h"
#include "algorithm/base.h"
#include "population.h"
#include "mpi_environment.h"

namespace pagmo
{

bool mpi_environment::m_initialised = false;
bool mpi_environment::m_multithread = false;

/// Default constructor.
/**
 * Initialises the MPI environment with MPI_Init_thread. pagmo::mpi_environment objects should be created only in the main
 * thread of execution.
 * 
 * @throws std::runtime_error if another instance of this class has already been created,
 * or if the MPI implementation does not support at least the MPI_THREAD_SERIALIZED thread level and this is the root node,
 * or if the world size is not at least 2.
 */
mpi_environment::mpi_environment()
{
	if (m_initialised) {
		pagmo_throw(std::runtime_error,"cannot re-initialise the MPI environment");
	}
	m_initialised = true;
	int thread_level_provided;
	MPI_Init_thread(NULL,NULL,MPI_THREAD_MULTIPLE,&thread_level_provided);
	if (thread_level_provided >= MPI_THREAD_MULTIPLE) {
		m_multithread = true;
	}
	if (get_rank()) {
		// If this is a slave, it will have to stop here, listen for jobs, execute them, and exit()
		// when signalled to do so.
		listen();
	}
	// If this is the root node, it will need to be able to call MPI from multiple threads.
	if (thread_level_provided < MPI_THREAD_SERIALIZED && get_rank() == 0) {
		pagmo_throw(std::runtime_error,"the master node must support at least the MPI_THREAD_SERIALIZED thread level");
	}
	// World sizes less than 2 are not allowed.
	if (get_size() < 2) {
		pagmo_throw(std::runtime_error,"the size of the MPI world must be at least 2");
	}
}

/// Destructor.
/**
 * Will send a shutdown signal to all processes with nonzero rank and call MPI_Finalize().
 */
mpi_environment::~mpi_environment()
{
	// In theory this should never be called by the slaves.
	pagmo_assert(!get_rank());
	pagmo_assert(m_initialised);
	std::pair<boost::shared_ptr<population>,algorithm::base_ptr> shutdown_payload;
	for (int i = 1; i < get_size(); ++i) {
		// Send the shutdown signal to all slaves.
		send(shutdown_payload,i);
	}
	MPI_Finalize();
	m_initialised = false;
}

void mpi_environment::check_init()
{
	if (!m_initialised) {
		pagmo_throw(std::runtime_error,"MPI environment has not been initialised");
	}
}

/// Probe for message.
/**
 * This method is thread-safe only if mpi_environment::is_multithread returns true.
 * 
 * @param[in] source rank of the processor that will be probed.
 * 
 * @return true if source sent a message to this process, false otherwise.
 * 
 * @throws std::runtime_error if the MPI environment has not been initialised.
 */
bool mpi_environment::iprobe(int source)
{
	check_init();
	MPI_Status status;
	int flag;
	MPI_Iprobe(source,0,MPI_COMM_WORLD,&flag,&status);
	return flag;
}

/// MPI world size.
/**
 * This method is thread-safe only if mpi_environment::is_multithread returns true.
 * 
 * @return the MPI world size.
 * 
 * @throws std::runtime_error if the MPI environment has not been initialised.
 */
int mpi_environment::get_size()
{
	check_init();
	int retval;
	MPI_Comm_size(MPI_COMM_WORLD,&retval);
	return retval;
}

/// MPI rank.
/**
 * This method is thread-safe only if mpi_environment::is_multithread returns true.
 * 
 * @return the MPI rank of the process.
 * 
 * @throws std::runtime_error if the MPI environment has not been initialised.
 */
int mpi_environment::get_rank()
{
	check_init();
	int retval;
	MPI_Comm_rank(MPI_COMM_WORLD,&retval);
	return retval;
}

/// Thread-safety of the MPI implementation.
/**
 * This method is always thread-safe.
 * 
 * @return true if the MPI implementation of this process is completely thred-safe (i.e., the
 * MPI implementation supports MPI_THREAD_MULTIPLE), false for any other MPI thread level.
 * 
 * @throws std::runtime_error if the MPI environment has not been initialised.
 */
bool mpi_environment::is_multithread()
{
	check_init();
	return m_multithread;
}

void mpi_environment::listen()
{
	std::pair<boost::shared_ptr<population>,algorithm::base_ptr> payload;
	while (true) {
		// Receive the payload from the master.
		recv(payload,0);
		// If the payload is empty, it is the shutdown payload.
		if (payload.first.get() == 0) {
			pagmo_assert(payload.second.get() == 0);
			break;
		}

		// Create copy of data received
		const boost::shared_ptr<population> pop_copy(new population(*payload.first));
		try {
			// Perform the evolution.
			payload.second->evolve(*payload.first);
		} catch (const std::exception &e) {
			std::cout << "MPI Remote Error during island evolution using " << payload.second->get_name() << ": " << e.what() << std::endl;
		} catch (...) {
			std::cout << "MPI Remote Error during island evolution using " << payload.second->get_name() << ", unknown exception caught. :(" << std::endl;
		}

		try {                
			// Send back to the master the evolved population.
			send(payload.first,0);
		} catch (const boost::archive::archive_exception &e) {
			std::cout << "MPI Send Error during island evolution using " << payload.second->get_name() << ": " << e.what() << std::endl;
			// Send back to the master the original population.
			send(pop_copy,0);
		} catch (...) {
			std::cout << "MPI Send Error during island evolution using " << payload.second->get_name() << ", unknown exception caught. :(" << std::endl;
			// Send back to the master the original population.
			send(pop_copy,0);
                }
	}
	// Destroy the MPI environment before exiting.
	MPI_Finalize();
	std::exit(0);
}

}
