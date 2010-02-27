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

// 16/05/2008: Initial version by Dario Izzo.

#ifndef PAGMO_ALGORITHM_PSO_H
#define PAGMO_ALGORITHM_PSO_H

#include <iostream>
#include <string>

#include "../../config.h"
#include "../basic/population.h"
#include "base.h"

namespace pagmo
{
namespace algorithm {

/// Particle Swarm Optimisation Solver (PSO)
/**
 * The canonic Particle Swarm Optimisation algorithm. Pseudo algorithm taken from
 * http://en.wikipedia.org/wiki/Particle_swarm_optimization [July 2008] and implemented for the class algorithm::base
 */

class __PAGMO_VISIBLE pso: public base
{
	public:
		/// Constructor (deprecated)
		/**
		 * Instantiates a PSO algorithm with strategy 3 (for backcompatibility with older pagmos)
		 * \param[in] gen Generation to be evolved
		 * \param[in] inertia The particle inertia
		 * \param[in] cognitive The particle cognitive component
		 * \param[in] social The particle social component
		 * \param[in] vcoeff Defines the initial particle velocity. Must be in [0,1]. When 0 the particle initial velocity is
		 * zero, when one it is a random vector between the lower and upper bounds
		*/
		pso(int gen, const double &inertia, const double &cognitive, const double &social, const double &vcoeff);

		/// Constructor (deprecated)
		/**
		 * Instantiates a PSO algorithm with default parameters and strategy 3 (for backcompatibility with older pagmos)
		 * \param[in] gen Generation to be evolved
		*/
		pso(int gen);

		/// Constructor
		/**
		 * Instantiates a PSO algorithm giving full control over the parameters
		 * \param[in] gen Generation to be evolved
		 * \param[in] inertia The particle inertia
		 * \param[in] cognitive The particle cognitive component
		 * \param[in] social The particle social component
		 * \param[in] vcoeff Defines the initial particle velocity. Must be in [0,1]. When 0 the particle initial velocity is
		 * zero, when one it is a random vector between the lower and upper bounds
		 * \param[in] strategy Defines the PSO startegy to be used to update the velocities
		*/
		pso(int gen, const double &inertia, const double &cognitive, const double &social, const double &vcoeff, const int &strategy);

		/// Constructor
		/**
		 * Instantiates a PSO algorithm with default parameters
		 * \param[in] gen Generation to be evolved
		 * \param[in] strategy Defines the PSO startegy to be used to update the velocities
		*/
		pso(int gen, const int &strategy);

		/// Algorithm
		/**
		 * It performs a call to the PSO algorithm evolving the population for gen generations
		 * \param[in] popin Starting population
		 * \return Evolved population
		*/
		virtual population evolve(const population &popin) const;
		virtual pso *clone() const {
			return new pso(*this);
		}

		virtual std::string id_object() const;
	private:
		virtual void log(std::ostream &) const;
		size_t generations;
		double omega;
		double eta1;
		double eta2;
		double vcoeff;
		int strategy;
};

}
}

#endif
