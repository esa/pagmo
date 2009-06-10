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

// 03/02/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_NM_ALGORITHM_H
#define PAGMO_NM_ALGORITHM_H

#include <iostream>
#include <utility>
#include <vector>

#include "../../../config.h"
#include "../basic/population.h"
#include "../problems/GOproblem.h"
#include "go_algorithm.h"

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
class __PAGMO_VISIBLE nm_algorithm: public go_algorithm {
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
		nm_algorithm(int, const double &, const double &, const double &, const double &);
		nm_algorithm(int);
		virtual Population evolve(const Population &) const;
		/// Clone method.
		virtual nm_algorithm *clone() const {return new nm_algorithm(*this);}
		/// Return unique string identifying the algorithm.
		virtual std::string id_object() const {return id_name(); }
	private:
		std::vector<double> center_mass(const simplex &) const;
		std::vector<double> sub_mult_add(const vertex &, const vertex &,
			const double &, const vertex &) const;
		void check_bounds(vertex &, const GOProblem &) const;
		double simplex_diameter(const simplex &) const;
		void shuffle_simplex(simplex &, const GOProblem &) const;
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

#endif
