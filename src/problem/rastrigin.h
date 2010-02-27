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

#ifndef PAGMO_PROBLEM_RASTRIGIN_H
#define PAGMO_PROBLEM_RASTRIGIN_H

#include "../config.h"
#include "../types.h"

namespace pagmo{ namespace problem {

/// The Rastrigin problem.
/**
 * \image html rastrigin.png "Two-dimensional Rastrigin function."
 * \image latex rastrigin.png "Two-dimensional Rastrigin function." width=5cm
 *
 * The objective function here is the generalised n-dimensional Rastrigin function:
 * \f[
 * 	F\left(x_1,\ldots,x_n\right) = 10 \cdot n + \sum_{i=1}^n x_i^2 - 10\cdot\cos\left( 2\pi \cdot x_i \right), \quad x_i \in \left[ -5.12,5.12 \right].
 * \f]
 * The global minimum is in the origin, where \f$ F\left( 0,\ldots,0 \right) = 0 \f$.
 *
 * @author Dario Izzo (dario.izzo@esa.int)
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE rastrigin: public base
{
	public:
		rastrigin(int);
		base_ptr clone() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
};

}}

#endif
