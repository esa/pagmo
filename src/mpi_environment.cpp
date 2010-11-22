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
#include <boost/mpi/environment.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>
#include <cstdlib>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "exceptions.h"
#include "algorithm/base.h"
#include "population.h"
#include "mpi_environment.h"

namespace pagmo
{

static int mpi_int_init = 0;
static char **mpi_char_init = 0;

mpi_environment::mpi_environment():m_environment(new boost::mpi::environment(mpi_int_init,mpi_char_init))
{
	boost::mpi::communicator world;
	if (world.rank()) {
		// If this is a slave, it will have to stop here, listen for jobs, execute them, and exit()
		// when signalled to do so.
// std::cout << "hello I'm slave rank " << world.rank() << '\n';
		listen();
	}
	// If this is the root node, just finish the construction.
// std::cout << "hello I'm the master\n";
}

mpi_environment::~mpi_environment()
{
	int rank, numtasks;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
// std::cout << "shutting down\n";
	// In theory this should never be called by the slaves.
	pagmo_assert(!rank);
	for (int i = 1; i < numtasks; ++i) {
		// Send the shutdown signal to all slaves.
		MPI_Send(0,0,MPI_CHAR,i,2,MPI_COMM_WORLD);
	}
// std::cout << "shut down\n";
}

int mpi_environment::size() const
{
	return boost::mpi::communicator().size();
}

void mpi_environment::listen()
{
	MPI_Status status;
	int flag, rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	while (true) {
		// Query and catch, if any, the shutdown message.
		MPI_Iprobe(0,2,MPI_COMM_WORLD,&flag,&status);
		if (flag) {
// std::cout << "seen shutdown signal " << rank << '\n';
			MPI_Recv(0,0,MPI_CHAR,0,2,MPI_COMM_WORLD,&status);
// std::cout << "received shutdown signal " << rank << '\n';
			break;
		}
		// Query and catch, if any, the payload to be evolved.
		MPI_Iprobe(0,0,MPI_COMM_WORLD,&flag,&status);
		if (flag) {
			// First receive the size.
			int size;
// std::cout << "slave receiving size " << rank << '\n';
			MPI_Recv(static_cast<void *>(&size),1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
// std::cout << "slave received size " << rank << '\n';
			// Prepare the vector of chars.
			std::vector<char> buffer_char1(boost::numeric_cast<std::vector<char>::size_type>(size),0);
// std::cout << "slave receiving payload " << rank << '\n';
			// Receive the payload.
			MPI_Recv(static_cast<void *>(&buffer_char1[0]),size,MPI_CHAR,0,1,MPI_COMM_WORLD,&status);
// std::cout << "slave received payload " << rank << '\n';
			// Build the string from the vector.
			const std::string buffer_str1(buffer_char1.begin(),buffer_char1.end());
			// Unpickle the payload.
			std::pair<boost::shared_ptr<population>,algorithm::base_ptr> in;
			std::stringstream ss1(buffer_str1);
			boost::archive::text_iarchive ia(ss1);
			ia >> in;
			// Actually perform the evolution.
			in.second->evolve(*in.first);
			std::stringstream ss2;
			boost::archive::text_oarchive oa(ss2);
			const boost::shared_ptr<population> out_pop(new population(*in.first));
			const algorithm::base_ptr out_algo = in.second->clone();
			const std::pair<const boost::shared_ptr<population>, const algorithm::base_ptr> out(out_pop,out_algo);
			oa << out;
			const std::string buffer_str2(ss2.str());
			std::vector<char> buffer_char2(buffer_str2.begin(),buffer_str2.end());
			// Send the size.
			size = boost::numeric_cast<int>(buffer_char2.size());
// std::cout << "slave sending size " << rank << '\n';
			MPI_Send(static_cast<void *>(&size),1,MPI_INT,0,0,MPI_COMM_WORLD);
// std::cout << "slave sent size " << rank << '\n';
			// Send the string.
// std::cout << "slave sending payload " << rank << '\n';
			MPI_Send(static_cast<void *>(&buffer_char2[0]),size,MPI_CHAR,0,1,MPI_COMM_WORLD);
// std::cout << "slave sent payload " << rank << '\n';
		} else {
			// Sleep a bit if there is nothing to do.
			boost::this_thread::sleep(boost::posix_time::milliseconds(10));
		}
	}
	// Destroy the MPI environment before exiting.
	m_environment.reset(0);
// std::cout << "alright, exiting\n";
	std::exit(0);
}

}
