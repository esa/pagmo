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

// 18/05/2008: Initial version by Dario Izzo.

#ifndef PAGMO_ALGORITHM_DE_H
#define PAGMO_ALGORITHM_DE_H

#include <iostream>
#include <string>

#include "../../config.h"
#include "../basic/population.h"
#include "base.h"

namespace pagmo
{
namespace algorithm {

/// Differential Evolution Solver (DE)
/**
 * The Differential Evolution solver by Storn. Implementation is taken from
 * http://www.icsi.berkeley.edu/~storn/code.html and modified for the algorithm::base class. The following
 * description also is taken from that web pages: Differential Evolution grew out of Ken Price's attempts
 * to solve the Chebychev Polynomial fitting Problem that had been posed to him by Rainer Storn.
 * A breakthrough happened, when Ken came up with the idea of using vector differences for perturbing the vector
 * population. Since this seminal idea a lively discussion between Ken and Rainer and endless ruminations and computer
 * simulations on both parts yielded many substantial improvements which make DE the versatile and robust tool
 * it is today. The "DE community" has been growing since the early DE years of 1994 - 1996 and ever more
 * researchers are working on and with DE.
 */

class __PAGMO_VISIBLE de: public base
{
	public:
		/// Constructor
		/**
		 * It instantiate a DE algorithm.
		 * \param[in] gen Generation to be evolved
		 * \param[in] F Scaling parameter
		 * \param[in] CR Crossover probability
		*/
		de(int gen, const double &F, const double &CR, int strategy);

		/// Algorithm
		/**
		 * It performs a call to the DE algorithm evolving the population for gen generations
		 * \param[in] popin Starting population
		 * \return Evolved population
		*/
		virtual population evolve(const population &popin) const;
		virtual de *clone() const {
			return new de(*this);
		}
		virtual std::string id_object() const;
	private:
		virtual void log(std::ostream &) const;
		const size_t	generations;
		const double	F;
		const double	CR;
		const int	strategy;
};

}
}

#endif
