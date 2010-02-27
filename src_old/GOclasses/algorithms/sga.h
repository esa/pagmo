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

// Created by Dario Izzo on 10/05/08.

#ifndef PAGMO_ALGORITHM_SGA_H
#define PAGMO_ALGORITHM_SGA_H

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../../config.h"
#include "../basic/population.h"
#include "base.h"

namespace pagmo
{
namespace algorithm {

/// A simple genetic algorithm (SGA)
/**
 * This is a basic implementation of a genetic algorithm. Selection is either roulette wheel or 20% best repeated 5 times.
 * Mutation is either random, gaussian or bounded. Crossover is exponential. Elitism is implemented by inserting the best
 * every given generations.
 */
class __PAGMO_VISIBLE sga: public base
{
	public:
		/// Constructor.
		/**
		 * Creates the SGA algorithm with a random mutation strategy and a roulette wheel selection
		 * \param[in] gen Generation number
		 * \param[in] CR Crossover rate (long chromosomes need higer rates)
		 * \param[in] M Mutation rate
		* \param[in] best Elitism..... every best generation best-so-far is reinserted
		 */
		sga(int gen, const double &CR, const double &M, int best);
		/// Constructor.
		/**
		 * Creates the SGA algorithm with a random mutation strategy and a roulette wheel selection
		 * \param[in] gen Generation number
		 * \param[in] CR Crossover rate (long chromosomes need higer rates)
		 * \param[in] M Mutation rate
		* \param[in] best Elitism..... every best generation best-so-far is reinserted
		* \param[in] mutationRange If bounded mutation is selected it regulates the range within which mutation occur
		* \param[in] mutationType (0-gaussian, 1-bounded, 2-random)
		* \param[in] selectionType (0-20% best, 1-roulette)
		 */
		sga(int gen, const double &CR, const double &M, int best, double mutationRange, int mutationType, int selectionType);
		virtual population evolve(const population &) const;
		virtual sga *clone() const {
			return new sga(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual void log(std::ostream &) const;
		const size_t	generations;

		const double 	CR;		//crossover
		const double	M;		//mutation
		const int		insert_best;
		double		MR;    	//mutation range
		int		MType;
		//mutation type: 0: Original Gaussian,
		//		 1: Bounded mutation value,
		//		 2: Random mutation value.
		int		SType;
		//selection type: 0: Best %20 x 5
		//		  1: Roulette selection
		mutable			rng_uint32 rng;
};

}
}

#endif
