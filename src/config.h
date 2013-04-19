/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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
PaGMO offers a generalization of the island model paradigm working for global and local optimization
algorithms. Its main parallelization approach makes use of multiple threads, but MPI is also implemented and can be 
mixed in with multithreading. It provides a set of C++ classes and their exposition in Python language (PyGMO)
as to allow the user to solve, in a parallel fashion, global optimization tasks in the form: \n\n
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

.... yes it is that good ... its framework is applicable to single-objective, multiple-objectives, 
continuous, integer, box-constrained, non linear constrained, stochastic, deterministic optimization!!!!

 * @see Izzo, D., PyGMO and PyKEP: Open Source Tools for Massively Parallel Optimization in Astrodynamics (the case of interplanetary trajectory optimization), International Conference on Astrodynamics Tools and Techniques - ICATT, 2012. 
 * @see Izzo, D., Rucinski, M., and Biscani, F., The Generalized Island Model, Parallel Architectures and Bioinspired Algorithms, Springer Berlin/Heidelberg, pp.151--169, 2012.
 * @see Rucinski, M., Izzo, D., and Biscani, F., On the Impact of the Migration Topology on the Island Model, Parallel Computing, 36, Elsevier, pp.555-571, 2010.
 
\section algorithms Algorithms currently implemented in PaGMO

A number of algorithms are already implemented in PaGMO and thus immediately
available after installation. The user can implement his own algorithms both in C++
and directly in Python. Check the algorithm documentation to verify whether a particolar problem class
can be fed to it (i.e. box-constrained, mixed integer etc.):

- A Simple Genetic Algorithm (pagmo::algorithm::sga)
- Differential Evolution [the original algorithm] (pagmo::algorithm::de)
- jDE [self-adaptive Differential volution] (pagmo::algorithm::jde)
- mde_pbx [self-adaptive Differential volution] (pagmo::algorithm::mde_pbx)
- Particle Swarm Optimization (pagmo::algorithm::pso)
- Particle Swarm Optimization [generational version] (pagmo::algorithm::pso_generational)
- Non-dominated Sorting Genetic Algorithm (pagmo::algorithm::nsga2)
- Adaptive Neighbourhood Simulated Annealing (pagmo::algorithm::sa_corana)
- Improved Harmony Search (pagmo::algorithm::ihs)
- Compass Search (pagmo::algorithm::cs)
- Monotonic Basin Hopping (pagmo::algorithm::mbh)
- Multistart (pagmo::algorithm::ms)
- Monte-Carlo (pagmo::algorithm::monte_carlo)

Other algorithm are available via third parties libraries, and can be included activating the respective
options in CMake, in particular:

- GSL library (open-source) -- includes Nelder-Mead, BFGS and more, see pagmo::algorithm::base_gsl
- NLOPT library (open-source) -- includes BOBYQA, COBYLA and more, see pagmo::algorithm::base_nlopt
- IPOPT library (open-source) -- includes IPOPT, see pagmo::algorithm::ipopt
- SNOPT library (commercial) -- includes SNOPT, see pagmo::algorithm::snopt

When working only in Python the SciPy algorithms from the 'optimize module' are available too (see http://docs.scipy.org/doc/scipy/reference/optimize.html).

\section problems Problems currently implemented in PaGMO

A number of global otimization problems are already implemented in PaGMO and thus
immediately available after installation. The user can implement his own problems both in C++ or
directly in Python.

- Classical Test Problems
 - Box-Constrained Contnuous Single-Objective
  - Paraboloid, Ackley, Rastrigin, Rosenbrock, Branin, Schwefel, Griewank, Lennard-Jones, Levy5, Himmelblau ....
 - Continuous, constrained Single-Objective
  - From Luksan-Vlcek book (3 problems also in the original IPOPT dist.), Toy-problem from SNOPT manual
 - Integer, constrained Single-Objective
  - Golomb Ruler, Knapsack Problem
 - Box-Constrained Integer Single-Objective
  - The string match problem
 - Multi-objective, continuous, box-constrained
  - The SCH, FON, DTLZ and ZDT problems 
 - Stochastic, continuous and box-constrained
  - A simple inventory problem
- Engineering Problems
  - All problems from the GTOP database (http://www.esa.int/gsp/ACT/inf/op/globopt.htm) and An Interplanetary, Multiple Gravity Assist, Low-Thrust problem (MGA-LT)
  - Spheres (pagmo::problem::spheres, pagmo::problem::spheres_q). Two problems of neurocontroller evolution for the MIT SPHERE test-bed on board the ISS
\section install Installation guide

To install PaGMO from source code you will need git and CMake installed in your system. 

On Unix systems:

- Clone the PaGMO git repository on your local machine: \code git clone git://git.code.sf.net/p/pagmo/code pagmo-code \endcode
- Create a build directory in your pagmo directory and move there: \code 
cd pagmo-code \endcode \code mkdir build \endcode \code cd build 
\endcode
- Run ccmake to configure your makefile (or project): \code ccmake ../ \endcode
- In ccmake, press c to configure, then (see figure below) select the options that are desired 
(e.g. compile the main file?, compile  PyGMO?) press c to configure again and then g to generate the makefile. 
Selecting the option PyGMO you will also build the python version of the code. 
In this case make sure you have python installed. CMake will try to locate the current installation 
directory of your python and install there the code upon make install.

\image html ccmake.png
- Build PaGMO: \code make \endcode
- Test PaGMO (if tests are enabled in ccmake): \code make test\endcode
- Install PaGMO: \code make install \endcode

On Windows systems, the procedure is analogous (you will likely use the Windows CMake GUI instead of ccmake). We have so 
far tested compilation only using MinGW

\section PyGMO Interactive python session

Our suggestion in using PaGMO is to activate the option PyGMO and, after installation, start an
interactive session with ipython. In the image below you see an example on how to solve the 100
dimensional pagmo::problem::schwefel using a pagmo::topology::ring evolution of 8 island (20 individuals each)
using pagmo::algorithm::de (Differential Evolution) in a 4 CPU machine.

\image html ipython.png

\section mpi_intro MPI support

By default PaGMO parallelizes the optimization process by opening multiple local threads of execution, and hence the parallelism is confined to a single machine. For use in cluster
environments, PaGMO can employ MPI (Message Passing Interface) to distribute the workload among multiple machines. Detailed instructions on how to enable and use the MPI support in PaGMO
can be found in \ref mpi_support "this page".
*/

/*!
\page serialization Serialization
blah
*/

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 100000 \
	+ __GNUC_MINOR__ * 1000 \
	+ __GNUC_PATCHLEVEL__ * 10)
#endif

#ifdef _WIN32
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
