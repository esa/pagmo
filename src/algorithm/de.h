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

#ifndef PAGMO_ALGORITHM_DE_H
#define PAGMO_ALGORITHM_DE_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "base.h"


namespace pagmo { namespace algorithm {

/// Differential Evolution Algorithm
/**
 *
 * \image html de.jpg "Differential Evolution block diagram."
 * \image latex de.jpg "Differential Evolution block diagram." width=5cm
 *
 * Differential Evolution is an heuristic optimizer developed by Rainer Storn and Kenneth Price.
 *
 * ''A breakthrough happened, when Ken came up with the idea of using vector differences for perturbing
 * the vector population. Since this seminal idea a lively discussion between Ken and Rainer and endless
 * ruminations and computer simulations on both parts yielded many substantial improvements which
 * make DE the versatile and robust tool it is today'' (from the official web pages....)
 *
 * The implementation provided for PaGMO derives from the code provided in the official
 * DE web site and is suitable for box-constrained single-objective continuous optimization.
 *
 * At each call of the evolve method a number of function evaluations equal to m_gen * pop.size()
 * is performed.
 *
 * NOTE: when called on mixed-integer problems de treats the integer part as fixed and optimizes
 * the continuous part.
 *
 * @see http://www.icsi.berkeley.edu/~storn/code.html for the official DE web site
 * @see http://www.springerlink.com/content/x555692233083677/ for the paper that introduces Differential Evolution
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */
		
class __PAGMO_VISIBLE de: public base
{
public:
	de(int, const double & = 0.8, const double & = 0.9, int = 2);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	// Number of generations.
	const int m_gen;
	// Weighting factor
	const double m_f;
	// Crossover probability
	const double m_cr;
	// Startegy
	const int m_strategy;
};

}} //namespaces

#endif // DE_H
