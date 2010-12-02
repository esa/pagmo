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

#include <boost/shared_ptr.hpp>
#include <cstdlib>
#include <mpi.h>
#include <utility>

#include "exceptions.h"
#include "algorithm/base.h"
#include "population.h"
#include "mpi_environment.h"

namespace pagmo
{

mpi_environment::mpi_environment()
{
	int thread_level_provided;
	MPI_Init_thread(NULL,NULL,MPI_THREAD_MULTIPLE,&thread_level_provided);
	if (get_rank()) {
		// If this is a slave, it will have to stop here, listen for jobs, execute them, and exit()
		// when signalled to do so.
		listen();
	}
	// If this is the root node, just finish the construction.
}

mpi_environment::~mpi_environment()
{
	// In theory this should never be called by the slaves.
	pagmo_assert(!get_rank());
	std::pair<boost::shared_ptr<population>,algorithm::base_ptr> shutdown_payload;
	for (int i = 1; i < get_size(); ++i) {
		// Send the shutdown signal to all slaves.
		send(shutdown_payload,i);
	}
	MPI_Finalize();
}

bool mpi_environment::iprobe(int source)
{
	MPI_Status status;
	int flag;
	MPI_Iprobe(source,0,MPI_COMM_WORLD,&flag,&status);
	return flag;
}

int mpi_environment::get_size()
{
	int retval;
	MPI_Comm_size(MPI_COMM_WORLD,&retval);
	return retval;
}

int mpi_environment::get_rank()
{
	int retval;
	MPI_Comm_rank(MPI_COMM_WORLD,&retval);
	return retval;
}

bool mpi_environment::is_multithread()
{
	int thread_level_provided;
	MPI_Query_thread(&thread_level_provided);
	return (thread_level_provided >= MPI_THREAD_MULTIPLE);
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
		// Perform the evolution.
		payload.second->evolve(*payload.first);
		// Send back to the master the evolved payload.
		send(payload,0);
	}
	// Destroy the MPI environment before exiting.
	MPI_Finalize();
	std::exit(0);
}

}
