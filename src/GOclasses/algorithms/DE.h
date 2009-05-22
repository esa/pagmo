/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

#ifndef PAGMO_DE_H
#define PAGMO_DE_H

#include <iostream>

#include "../../config.h"
#include "../basic/population.h"
#include "go_algorithm.h"

/// Differential Evolution Solver
/**
 * This implementation is taken from Storn original paper. Implementation is taken from
 * http://www.icsi.berkeley.edu/~storn/code.html
 */

class __PAGMO_VISIBLE DEalgorithm: public go_algorithm {
	public:
                /// Constructor
                /**
                 * It instantiate a DE algorithm.
                 * \param[in] gen Generation to be evolved
                 * \param[in] F DE scaling parameter (algorithm specific)
                 * \param[in] CR crossover probability (algorithm specific)
                */
                DEalgorithm(int gen, const double &F, const double &CR, int strategy);

                /// Algorithm
                /**
                 * It performs a call to the DE algorithm evolving the population for gen generations
                 * \param[in] popin Starting population
                 * \return Evolved population
                */
                virtual Population evolve(const Population &popin) const;
		virtual DEalgorithm *clone() const {return new DEalgorithm(*this);}
		virtual std::string id_object() const;
	private:
		virtual void log(std::ostream &) const;
		const size_t	generations;
		const double	F;
		const double	CR;
		const int		strategy;
};

#endif
