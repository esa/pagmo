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
#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>
#include <cstdlib>
#include <sstream>
#include <string>
#include <utility>

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
std::cout << "hello I'm slave rank " << world.rank() << '\n';
		listen();
	}
	// If this is the root node, just finish the construction.
std::cout << "hello I'm the master\n";
}

mpi_environment::~mpi_environment()
{
std::cout << "shutting down\n";
	boost::mpi::communicator world;
	// In theory this should never be called by the slaves.
	pagmo_assert(!world.rank());
	for (int i = 1; i < world.size(); ++i) {
		// Send the shutdown signal to all slaves.
		world.send(i,1);
	}
std::cout << "shut down\n";
}

void mpi_environment::listen()
{
	boost::mpi::communicator world;
	std::pair<boost::shared_ptr<population>,algorithm::base_ptr> in;
	while (true) {
		// Query and catch, if any, the shutdown message.
		if (world.iprobe(0,1)) {
std::cout << "seen shutdown signal " << world.rank() << '\n';
			world.recv(0,1);
std::cout << "received shutdown signal " << world.rank() << '\n';
			break;
		}
		// Query and catch, if any, the payload to be evolved.
		if (world.iprobe(0,0)) {
			std::string payload;
std::cout << "slave receiving " << world.rank() << '\n';
			world.recv(0,0,payload);
			std::stringstream ss1(payload);
			boost::archive::text_iarchive ia(ss1);
			ia >> in;
std::cout << "slave received " << world.rank() << '\n';
			in.second->evolve(*in.first);
std::cout << "slave sending " << world.rank() << '\n';
			std::stringstream ss2;
			boost::archive::text_oarchive oa(ss2);
			const boost::shared_ptr<population> out_pop(new population(*in.first));
			const algorithm::base_ptr out_algo = in.second->clone();
			const std::pair<const boost::shared_ptr<population>, const algorithm::base_ptr> out(out_pop,out_algo);
			oa << out;
			payload = ss2.str();
			world.send(0,0,payload);
std::cout << "slave sent " << world.rank() << '\n';
		}
		// Sleep a bit if there is nothing to do.
		boost::this_thread::sleep(boost::posix_time::milliseconds(100));
	}
	// Destroy the MPI environment before exiting.
	m_environment.reset(0);
std::cout << "alright, exiting\n";
	std::exit(0);
}

}
