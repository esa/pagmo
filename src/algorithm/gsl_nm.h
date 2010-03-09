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

#ifndef PAGMO_ALGORITHM_GSL_NM_H
#define PAGMO_ALGORITHM_GSL_NM_H

#include <cstddef>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "../config.h"
#include "../population.h"
#include "base.h"
#include "base_gsl.h"

namespace pagmo { namespace algorithm {

/// GSL Nelder-Mead wrapper.
/**
 * Wrapper around the implementation of the Nelder-Mead simplex method available in the GNU Scientific Library (GSL).
 * The GSL function used is called gsl_multimin_fminimizer_nmsimplex2.
 *
 * <b>Usage notes</b>: to increase the convergence of this algorithm, try increasing the number of maximum iterations performed.
 *
 * @see http://www.gnu.org/software/gsl/manual/html_node/Multimin-Algorithms-without-Derivatives.html
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE gsl_nm: public base, base_gsl
{
	public:
		gsl_nm(int max_iter = 100, const double &tol = 1E-6, const double &step_size = 1);
		base_ptr clone() const;
		void evolve(population &) const;
	protected:
		std::string human_readable_extra() const;
	private:
		static void cleanup(gsl_vector *, gsl_vector *, gsl_multimin_fminimizer *);
		static void check_allocs(gsl_vector *, gsl_vector *, gsl_multimin_fminimizer *);
	private:
		std::size_t	m_max_iter;
		double		m_tol;
		double		m_step_size;
};

}}

#endif
