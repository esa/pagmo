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

/*! \mainpage PaGMO (Parallel Global Multiobjective Optimizer)
\section intro_sec Introduction 
PaGMO is an implementation of the island model paradigm to parallelize global and local optimization algorithms using multiple threads (tested to up to 10000 threads). It provides a set of C++ classes
and their exposition in Python language as to allow the user to solve, in a parallel fashion, global optimization tasks in the form: \n\n
<center>
  \f$
    \begin{array}{rl}
    \mbox{find:} & \mathbf x \in R^n \\
    \mbox{to minimize:} & \mathbf f(\mathbf x, \mathbf s) \\
    \mbox{subject to:} & \mathbf {lb} \le \mathbf x \le \mathbf {ub}
    \end{array}
  \f$
</center>\n\n

\section algorithms Algorithms currently implemented in PaGMO

A number of algorithms are already implemented in PaGMO and thus immediately available after installation. The user can implement his own algorithms only in C++ and recompiling the code
 
- A simple genetic algorithm (SGA)
- Differential Evolution (DE)
- Particle Swarm Optimization (PSO)
- Multiple Particle Swarm optimization (MPSO)
- Adaptive Neighbourhood Simulated Annealing (AN-SA)
- Improved Harmony Search (IHS)
- Nelder-Mead (NM)
- Compass Search (CS)

In the next future, the inclusion of the following algorithms is being discussed: SNOPT, IPOPT, Basin-Hopping, Simulated Annealing.
 
\section problems Problems currently implemented in PaGMO

A number of GO problems are already implemented in PaGMO and thus immediately available after installation. The user can implement his own problems only in C++ and recompiling the code.

- Classical Test Problems
 - Continuous and box bounded
  - Ackley
  - Rastrigin
  - Rosenbrock
  - Griewank
  - Lennard-Jones
  - Levy
 - Stochastic, continuous and box bounded
  - Inventory Problem
- Engineering Problems
  - Interplanetary Trajectories (from the GTOP database)
   - MGA type of problems
    - GTOC1 
    - Cassini1
   - MGA-DSM type of problems
    - Messenger, MessengerFull
    - Cassini2
    - TandEM
    - Laplace
    - Rosetta
   - Low-Thrust type of problems
    - Earth-Mars (Impulsive Transcription)
    - Earth-Mars (Constant-Thrust Transcription)
  - Evolutionary Robotics Problems
    - Twodee (requires thrid-party console application)
    - Marsrover (requires third-party console application)


\section install Installation guide

To install PaGMO from source code you will need git and cmake, ccmake installed in your system.

- Clone the PaGMO git repository on your local machine: \code git clone git://pagmo.git.sourceforge.net/gitroot/pagmo/pagmo \endcode
- Create a build directory in your pagmo directory and move there: \code cd pagmo \endcode \code mkdir build \endcode \code cd build \endcode
- Run cmake to configure your makefile (or project): \code ccmake ../ \endcode
- In ccmake, press c to configure, then (see figure below) select the options that are desired (e.g. compile the main file?, compile  PyGMO?) press c to configure again and then g to generate the makefile. Selecting the option PyGMO you will also build the python version of the code. In this case make sure you have python installed. Cmake will try to locate the current installation directory of your python and install there the code.

\image html ccmake.png
- Build PaGMO: \code make \endcode
- Install PaGMO: \code make install \endcode

\section PyGMO Interactive python session

Our suggestion in using PaGMO is to activate the option PyGMO and, after installation, start an interactive session with ipython. In the image below you see an example on how to solve the 100 dimensional ackley problem using a ring topology parallel evolution of 8 island using differential evolution

\image html ipython.png
*/
#ifndef PAGMO_CONSTANTS_H
#define PAGMO_CONSTANTS_H

#ifdef __GNUC__

#include <cmath>

#else

#define M_PI (3.1415926535897931)
#define M_PI_4 (0.78539816339744828)
#define M_PI_2 (1.5707963267948966)

#endif

#endif
