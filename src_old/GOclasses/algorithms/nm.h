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

// 03/02/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_ALGORITHM_NM_H
#define PAGMO_ALGORITHM_NM_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../../config.h"
#include "../basic/population.h"
#include "../problems/base.h"
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
 */
class __PAGMO_VISIBLE nm: public base
{
		/// Vertex type.
		/**
		 * A vertex is analogous to a decision vector.
		 */
		typedef std::vector<double> vertex;
		/// Simplex type.
		/**
		 * A simplex is defined here as a vector of vertex - fitness pairs. We explicitly
		 * store the fitness in order to avoid objective function evaluations.
		 */
		typedef std::vector<std::pair<vertex,double> > simplex;
		struct sorter;
	public:
		nm(int, const double &, const double &, const double &, const double &);
		nm(int);
		virtual population evolve(const population &) const;
		/// Clone method.
		virtual nm *clone() const {
			return new nm(*this);
		}
		/// Return unique string identifying the algorithm.
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		std::vector<double> center_mass(const simplex &) const;
		std::vector<double> sub_mult_add(const vertex &, const vertex &,
		                                 const double &, const vertex &) const;
		void check_bounds(vertex &, const problem::base &) const;
		double simplex_diameter(const simplex &) const;
		void shuffle_simplex(simplex &, const problem::base &) const;
		virtual void log(std::ostream &) const;
		/// Number of generations.
		const size_t	m_gen;
		/// Reflection coefficient (\f$ \alpha \f$).
		const double	m_alpha;
		/// Expansion coefficient (\f$ \gamma \f$).
		const double	m_gamma;
		/// Contraction coefficient (\f$ \rho \f$).
		const double	m_rho;
		/// Shrink coefficient (\f$ \sigma \f$).
		const double	m_sigma;
};

}
}

#endif
