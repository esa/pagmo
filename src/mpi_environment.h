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

#ifndef PAGMO_MPI_ENVIRONMENT_H
#define PAGMO_MPI_ENVIRONMENT_H

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/utility.hpp>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

#include "config.h"

/*!
\page mpi_support MPI support in PaGMO
PaGMO can be configured to use MPI to dispatch parallel optimisation instances to computers participating in a cluster.

@author Francesco Biscani (bluescarni@gmail.com)
@author Dante Stroe (dante.stroe@gmail.com)

\section mpi_requirements Requirements
In order to enable and use MPI support in PaGMO a standard-compliant MPI implementation (e.g., Open MPI, MPICH2, etc.) must be available
on all systems participating to the cluster. PaGMO uses basic MPI 1.2 calls such as MPI_Recv, MPI_Send etc., and does not employ any function
specific to MPI 2.x.

For best performance, the root node of the MPI cluster (i.e., the node where mpiexec/mpirun is launched) should have a thread-safe
MPI implementation. This means that a call to MPI_Query_thread should ideally return MPI_THREAD_MULTIPLE. It is possible to run PaGMO
also if MPI_Query_thread returns MPI_THREAD_SERIALIZED, but in this case there will be a performance penalty. With any lesser level of thread
support in the root node, PaGMO will refuse to operate in MPI mode.

@see http://www.fz-juelich.de/jsc/juropa/www3/MPI_Query_thread.html

In a PaGMO MPI cluster, problem, population and algorithm instances are sent from the master node to the slave nodes where the optimisation process
takes place. In order for the objects to be sent over the network, they need to be serializable. Serialization for all classes shipped with PaGMO is
already implemented. The writer of new PaGMO classes will need to make sure that her own classes are serializable. 

Additionally, it must always kept in mind that when working with MPI each processor is living inside a separate process: aside from the objects serialized
and dispatched around in the cluster, there is no other form of communication. For instance, problems or algorithms employing a shared memory state
(e.g., static class members) must be carefully treated during serialization in order to provide to the remote processes all the information needed
for the optimisation process.
\section mpi_compilation Compilation
In order to compile PaGMO with MPI support, just enable the ENABLE_MPI option in CMake. The build system should automatically detect your MPI implementation
and configure the build accordingly. In case of issues, you can explictly specify the location of you MPI libraries and headers by manipulating the relevant
CMake variables.
\section mpi_model PaGMO's MPI model
MPI support in PaGMO is implemented via the pagmo::mpi_island and pagmo::mpi_environment classes. The pagmo::mpi_environment class is used just for initialisation purposes,
and an instance of it should be created at the very beginning of any MPI-enabled PaGMO main() function.

The pagmo::mpi_island class is an island class deriving from pagmo::base_island that can be operated exactly like any other island class (e.g., pagmo::island):
it can be created and operated stand-alone, or inserted in a pagmo::archipelago. Different types of islands can exist in the same pagmo::archipelago, so that it is
possible to mix MPI islands with islands of different kind.

At the beginning of a PaGMO MPI main() execution, a list of available processors is created. The list contains the integers from 1 to N - 1, where N is the size
of the MPI world (from which it follows that it is not possible to use PaGMO with an MPI world of size less than 2). In the PaGMO model, the process with rank 0
is tasked with the dispatch and coordination of jobs and does not perform any optimization process.

Whenever an evolution method is called from from a pagmo::mpi_island, the island will check the list of available processors and, if a processor is available, will
erase the processor ID from the list and dispatch the evolution to that processor. At the end of the evolution, the island retrieves the payload and adds the processor ID back
to the list of available processors.

Whenever the number of MPI islands is at least equal to the MPI world size, it might happen that one or more islands are not able to acquire any processor at the beginning of
the evolution, all the processors being busy. In such case a fair priority queue is created, and the islands waiting for a processor to be released are added to the
end of the queue. Whenever a processor is released, the queue is notified and the first island in the queue acquires the processor and procedes as above.
\section mpi_example MPI example
The following simple example shows how to create and use MPI islands in a PaGMO C++ main().
\code
// Include the global PaGMO header.
#include "pagmo.h"

using namespace pagmo;

int main()
{
	// Initialise the MPI environment.
	mpi_environment env;
	// Create a problem and an algorithm.
	problem::dejong prob(10);
	algorithm::monte_carlo algo(100);
	// Create an archipelago of 10 MPI islands.
	archipelago a;
	a.set_topology(topology::ring());
	for (int i = 0; i < 10; ++i) {
		a.push_back(mpi_island(algo,prob,1));
	}
	// Evolve the archipelago 10 times.
	a.evolve(10);
	a.join();
	return 0;
}
\endcode
The code is exactly the same as it would be for use with pagmo::island instances, the only difference being the creation of a pagmo::mpi_environment at the very beginning.
\section mpi_execution Executing MPI programs
An MPI-enabled PaGMO executable can be executed just like any MPI executable. The following instructions assume that
the MPI environment is from Open MPI; different MPI implementations might need slightly different setup (e.g., MPICH2
requires a daemon to be running for dispatching MPI jobs).
\subsection local_mpi_execution Local execution
Local MPI execution is a simple way to verify that the MPI program is working correctly. As the name implies,
all jobs will be started on the local machine.

Assuming the name of the executable is 'main', local MPI execution can be launched with:
\verbatim
$ mpiexec -n 10 ./main
\endverbatim
where the argument to the
\verbatim
-n
\endverbatim
parameter is the MPI world size (i.e., the number of MPI jobs that
will be launched concurrently - in this case 10).
\subsection cluster_mpi_execution Execution in a cluster
Once it has been verified that local execution works as expected, the next step is to run the MPI-enabled
PaGMO executable in a cluster. Again, the setup is the same needed for any MPI executable. Namely:

  - exactly the same executable should be present on all nodes of the cluster;
  - the MPI executable should reside in the same path on all nodes;
  - any resource needed by the executable (e.g., files) should be available at the same location
    on all nodes.

To run the executable on the cluster, a hostfile will be needed. A hostfile is a plain text file describing
the machine participating in the cluster. A sample sample hostfile is the following:
\verbatim
sophia.estec.esa.int slots=8
ursula.estec.esa.int slots=8
\endverbatim
For each line, the first entry is the name (or IP address) of a computer participating in the cluster, the second
entry is the number of jobs that can be run on that machine. Typically the slots number will be equal to the number
of CPUs on the computer.

Once the hostfile has been created, the MPI executable can be launched with the command:
\verbatim
$ mpiexec -n 10 -hostfile hfile.txt ./main
\endverbatim
where hfile.txt is your hostfile.
*/

namespace pagmo
{

/// MPI environment class.
/**
 * This class is used to initialise the PaGMO MPI environment: an instance of this class should be created
 * before using any MPI feature in PaGMO. Apart from this kind of usage, regular users should never need to
 * access any method from this class. See the \ref mpi_support "MPI page" for a usage example.
 * 
 * <b>NOTE</b>: this class is available only if PaGMO was compiled with MPI support.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE mpi_environment: private boost::noncopyable
{
	public:
		mpi_environment();
		~mpi_environment();
		static bool is_multithread();
		static int get_size();
		static int get_rank();
		/// Receive MPI payload.
		/**
		 * Receive an instance of class T from the processor with ID source and store it into retval.
		 * This method is thread-safe only if mpi_environment::is_multithread returns true.
		 * 
		 * @param[out] retval instance of class T that will contain the payload.
		 * @param[in] source rank of the processor from which the message will be received.
		 */
		template <class T>
		static void recv(T &retval, int source)
		{
			check_init();
			MPI_Status status;
			// First receive the size.
			int size;
			MPI_Recv(static_cast<void *>(&size),1,MPI_INT,source,0,MPI_COMM_WORLD,&status);
			// Prepare the vector of chars.
			std::vector<char> buffer_char(boost::numeric_cast<std::vector<char>::size_type>(size),0);
			// Receive the payload.
			MPI_Recv(static_cast<void *>(&buffer_char[0]),size,MPI_CHAR,source,1,MPI_COMM_WORLD,&status);
			// Build the string from the vector.
			const std::string buffer_str(buffer_char.begin(),buffer_char.end());
			// Unpickle the payload.
			std::stringstream ss(buffer_str);
			boost::archive::text_iarchive ia(ss);
			ia >> retval;
		}
		/// Send MPI payload.
		/**
		 * Send an instance of class T to the processor with ID destination.
		 * This method is thread-safe only if mpi_environment::is_multithread returns true.
		 * 
		 * @param[in] payload instance of class T that will be sent to destination.
		 * @param[in] destination rank of the processor to which the message will be sent.
		 */
		template <class T>
		static void send(const T &payload, int destination)
		{
			check_init();
			std::stringstream ss;
			boost::archive::text_oarchive oa(ss);
			oa << payload;
			const std::string buffer_str(ss.str());
			std::vector<char> buffer_char(buffer_str.begin(),buffer_str.end());
			// Send the size.
			int size = boost::numeric_cast<int>(buffer_char.size());
			MPI_Send(static_cast<void *>(&size),1,MPI_INT,destination,0,MPI_COMM_WORLD);
			// Send the string.
			MPI_Send(static_cast<void *>(&buffer_char[0]),size,MPI_CHAR,destination,1,MPI_COMM_WORLD);
		}
		static bool iprobe(int);
	private:
		static void listen();
		static void check_init();
		static bool	m_initialised;
		static bool	m_multithread;
};

}

#endif
