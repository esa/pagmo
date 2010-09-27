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
 *   the Free Software Foundation; either version 3 of the License, or       *
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

#ifndef PAGMO_CONFIG_H
#define PAGMO_CONFIG_H

/*! \mainpage PaGMO (Parallel Global Multiobjective Optimizer)
\section intro_sec Introduction 
PaGMO is a generalization of the island model paradigm to parallelize global and local optimization algorithms using multiple threads (tested to up to 10000 threads). It provides a set of C++ classes
and their exposition in Python language as to allow the user to solve, in a parallel fashion, global optimization tasks in the form: \n\n
<center>
	\f$
	\begin{array}{rl}
	\mbox{find:} & \mathbf x \in R^n \times N^m \\
	\mbox{to minimize:}	& \mathbf f(\mathbf x, \mathbf s) \\
	\mbox{subject to:}	& \mathbf {lb} \le \mathbf x \le \mathbf {ub} \\
				& \mathbf c(\mathbf x) = \mathbf 0 \\
				& \mathbf c_{in}(\mathbf x) \le \mathbf 0
	\end{array}
	\f$
</center>\n\n

\section algorithms Algorithms currently implemented in PaGMO

A number of algorithms are already implemented in PaGMO and thus immediately
available after installation. The user can implement his own algorithms both in C++
and directly in Python. Check the algorithm documentation to verify whether a particolar problem class
can be fed to it (i.e. box-constrained or mixed integer etc.):

- A Simple Genetic Algorithm (pagmo::algorithm::sga)
- Differential Evolution (pagmo::algorithm::de)
- Particle Swarm Optimization (pagmo::algorithm::pso)
- Adaptive Neighbourhood Simulated Annealing (pagmo::algorithm::sa_corana)
- Improved Harmony Search (pagmo::algorithm::ihs)
- Compass Search (pagmo::algorithm::cs)
- Monotonic Basin Hopping (pagmo::algorithm::mbh)
- Multistart (pagmo::algorithm::ms)
- Monte-Carlo (pagmo::algorithm::monte_carlo)

Other algorithm are available via third parties libraries, and can be included activating the respective
options in ccmake, in particular:

- GSL library (open-source) -- includes Nelder-Mead, BFGS and more, see pagmo::algorithm::base_gsl
- NLOPT library (open-source) -- includes bobyqa, cobyla and more, see pagmo::algorithm::base_nlopt
- IPOPT library (open-source) -- includes IPOPT, see pagmo::algorithm::ipopt
- SNOPT library (commercial) -- includes SNOPT, see pagmo::algorithm::snopt

When working only in Python the scipy algorithms are available too.

\section problems Problems currently implemented in PaGMO

A number of global otimization problems are already implemented in PaGMO and thus
immediately available after installation. The user can implement his own problems both in C++ or
directly in Python.

- Classical Test Problems
 - Continuous, box bounded
  - Paraboloid, Ackley, Rastrigin, Rosenbrock, Branin, Schwefel, Griewank, Lennard-Jones, Levy5, HimmelBlau
 - Continuous, constrained
  - From Luksan-Vlcek book (3 problems also in the original IPOPT dist.), Toy-problem from SNOPT manual
 - Integer Programming
  - Golomb Ruler, Knapsack Problem
 - Multi-objective, continuous, box-constrained
  - The SCH and FON problems from nsga-II
 - Stochastic, continuous and box bounded
  - xxx
- Engineering Problems
  - All problems from the GTOP database, An Interplanetary, Multiple Gravity Assist, Low-Thrust problem (MGA-LT)

\section install Installation guide

To install PaGMO from source code you will need git and cmake, ccmake installed in your system.

- Clone the PaGMO git repository on your local machine: \code git clone git://pagmo.git.sourceforge.net/gitroot/pagmo/pagmo \endcode
- Create a build directory in your pagmo directory and move there: \code cd pagmo \endcode \code mkdir build \endcode \code cd build \endcode
- Run cmake to configure your makefile (or project): \code ccmake ../ \endcode
- In ccmake, press c to configure, then (see figure below) select the options that are desired (e.g. compile the main file?, compile  PyGMO?) press c to configure again and then g to generate the makefile. Selecting the option PyGMO you will also build the python version of the code. In this case make sure you have python installed. Cmake will try to locate the current installation directory of your python and install there the code.

\image html ccmake.png
- Build PaGMO: \code make \endcode
- Test PaGMO (if tests are enabled in ccmake): \code make test\endcode
- Install PaGMO: \code make install \endcode

\section PyGMO Interactive python session

Our suggestion in using PaGMO is to activate the option PyGMO and, after installation, start an
interactive session with ipython. In the image below you see an example on how to solve the 100
dimensional pagmo::problem::schwefel using a pagmo::topology::ring evolution of 8 island (20 individuals each)
using pagmo::algorithm::de (Differential Evolution) in a 4 CPU machine.

\image html ipython.png

\section PyGMO MPI implementation 

The current MPI implementation works with C++ programs that make use of the pagmo libraries. In order to be able to run your solver on multiple processor you need to have a working version of MPI installed on your system. Any of the following MPI distribution should be compatibale with boost and implictly PAGMO:
OPENMPI, MPICH2, LAM/MPI. On debian systems the installation of these distributions is as simple as "sudo apt-get install mpich2"/"sudo apt-get install openmpi"
	Notes for solving approximation problems with PAGMO problems using MPI:
	The main difference is that in this case the archipelago are created using "mpi_islands" objects instead of "island", as follows:
@verbatim
a.push_back(pagmo::mpi_island(problem,algorithm,no_individuals,processor_id));
@endverbatim
	, where processor_id is an integer that represents the rank of the processor to which that particular island is assigned to.
	An archipelago with n mpi_islands can also be created by using the archipelago constructor, and setting the is_parallel attribute to "true" (default is "false"):
@verbatim
pagmo::archipelago b = pagmo::archipelago(problem, algorithm, n, n_indiv, topology, migration, direction, true)
@endverbatim
	, where "n" represents the number of islands of the archipelago, and the "is_parallel" is the boolean value that needs to be set to "true" in order to make use of the MPI environment. In case the number of islands is larger than the number of processors available, then the the computations of the islands' evoutions, are assigned to the processors in a round-robin manner.
	Note that the program must be treated as an MPI program (running on multiple processes/ors). While the initialization of the archipelago needs to be done on all processors, and the evolution of the islands called on all processors (the process filtering being done by the exposed "perform_evolution method" of the mpi_island class) other parts of the code like printing the results needs to be done onl  on the root process (the one with rank 0). An example can be found with the test_mpi.cpp program located in the "tests" folder.
	Testing: Once the pagmo base code, along with C++ MPI problem have been compiled, one can test their implementation by first starting the mpd deamon (on a single machine if it's tested with multiple processes on a single computer, or on all the machines that are being used). To run the problem one needs to execute
@verbatim
mpirun -np number_of_processes (-hostfile hostfilename)./mpi_program_executable
@endverbatim
, where the hostfile is used to specify the available machines in the case where the program is tested on multiple nodes instead of just one.
*/

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 100000 \
	+ __GNUC_MINOR__ * 1000 \
	+ __GNUC_PATCHLEVEL__ * 10)
#endif

#ifdef PAGMO_WIN32
	#ifdef PAGMO_DLL_EXPORT_API
		#define __PAGMO_VISIBLE __declspec(dllexport)
	#elif defined ( PAGMO_DLL_IMPORT_API )
		#define __PAGMO_VISIBLE __declspec(dllimport)
	#else
		#define __PAGMO_VISIBLE
	#endif
	#define __PAGMO_VISIBLE_FUNC __PAGMO_VISIBLE
#else
	#define __PAGMO_VISIBLE __attribute__ ((visibility("default")))
	#define __PAGMO_VISIBLE_FUNC
#endif

/// Root PaGMO namespace.
namespace pagmo {}

#endif
