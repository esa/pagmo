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
These pages contain the documentation of PaGMO (C++) API. PaGMO offers a generalization of the island model paradigm working for global and local optimization
algorithms. Its main parallelization approach makes use of multiple threads, but MPI is also implemented and can be 
mixed in with multithreading. PaGMO can be used to solve in a parallel fashion, global optimization tasks in the form: \n\n
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
 

\section mpi_intro MPI support

By default PaGMO parallelizes the optimization process by opening multiple local threads of execution, and hence the parallelism is confined to a single machine. For use in cluster
environments, PaGMO can employ MPI (Message Passing Interface) to distribute the workload among multiple machines. Detailed instructions on how to enable and use the MPI support in PaGMO can be found in \ref mpi_support "this page".
*/


// We solve the library madness in windows
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


// Redefining M_PI end M_E to avoid problems in windows
#include <boost/math/constants/constants.hpp>
#ifndef M_PI
#define M_PI boost::math::constants::pi<double>()
#endif
#ifndef M_PI_2
#define M_PI_2 boost::math::constants::pi<double>()/2
#endif
#ifndef M_PI_4
#define M_PI_4 boost::math::constants::pi<double>()/4
#endif
#ifndef M_E
#define M_E boost::math::constants::e<double>()
#endif

/// Root PaGMO namespace.
namespace pagmo {}

#endif
