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

#ifndef PAGMO_ALGORITHM_NELDER_MEAD_H
#define PAGMO_ALGORITHM_NELDER_MEAD_H

#include <cstddef>
#include <string>

#include "../config.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base.h"

namespace pagmo
{
namespace algorithm {

/// Nelder-Mead algorithm.
/**
 * \image html Nelder_Mead2.gif
 *
 * The Nelder-Mead method (or downhill simplex method) is a commonly used nonlinear
 * optimization algorithm by John Nelder and R. Mead (1965). It is a numerical method for minimizing an
 * objective function in a many-dimensional space.
 *
 * Each individual of the incoming population is seen as a vertex of an \f$ n+1 \f$-dimensional simplex (where \f$ n \f$ is the dimension
 * of the problem).
 *
 * The algorithm provided here implements the variant described in http://en.wikipedia.org/wiki/Nelder-Mead_method (as it appeared
 * on 27/05/2009).
 *
 * This algorithm is suitable for continuous, constrained and multiobjective optimisation, and will use the problem's ranking methods to establish whether
 * a decision vector is better than the other. Only the continuous part of the decision vectors will be optimised.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE nelder_mead: public base
{
	public:
		nelder_mead(int n_it, const double &thresh = 1E-10, const double &alpha = 1, const double &gamma = 2, const double &rho = .5, const double &sigma = .5);
		void evolve(population &) const;
		base_ptr clone() const;
	protected:
		std::string human_readable_extra() const;
	private:
		static void center_mass(decision_vector &, const population::const_iterator &, const population::const_iterator &,
			const problem::base::size_type &);
		static void sub_mult_add(decision_vector &, const decision_vector &, const decision_vector &, const double &,
			const decision_vector &, const problem::base::size_type &);
		static bool check_bounds(decision_vector &, const problem::base &);
		void simplex_size_check(const population::const_iterator &, const population::const_iterator &, population &, const problem::base::size_type &) const;
	private:
		// Number of generations.
		const std::size_t	m_gen;
		// Simplex rebuild threshold.
		const double		m_thresh;
		// Reflection coefficient.
		const double		m_alpha;
		// Expansion coefficient.
		const double		m_gamma;
		// Contraction coefficient.
		const double		m_rho;
		// Shrink coefficient.
		const double		m_sigma;
};

}
}

#endif
