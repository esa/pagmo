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

#ifndef PAGMO_ALGORITHM_GSL_GRADIENT_H
#define PAGMO_ALGORITHM_GSL_GRADIENT_H

#include <cstddef>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <string>

#include "../config.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base_gsl.h"

namespace pagmo { namespace algorithm {

/// Wrapper for GSL minimisers with derivatives
/**
 * This class can be used to build easily a wrapper around a GSL minimiser with derivatives. The gradient of the
 * objective function will be calculated numerically via the gsl_deriv_central GSL function. The minimisation is performed
 * by invoking the evolve_gradient() method.
 *
 * This class of minimisers supports single-objective, unconstrained, continuous optimisation.
 */
class __PAGMO_VISIBLE gsl_gradient: public base_gsl
{
	public:
		gsl_gradient(int, const double &, const double &, const double &, const double &);
	protected:
		void evolve_gradient(population &, const gsl_multimin_fdfminimizer_type *) const;
		std::string hr_extra() const;
	private:
		// Structure to feed parameters to the numerical differentiation wrapper.
		struct objfun_numdiff_wrapper_params
		{
			// Pointer to the problem.
			problem::base const		*prob;
			// Decision vector.
			decision_vector			x;
			// Fitness vector.
			fitness_vector			f;
			// Coordinate of the gradient being computed.
			problem::base::size_type 	coord;
		};
		static double objfun_numdiff_wrapper(double, void *);
		static void objfun_numdiff_central(gsl_vector *, const problem::base &, const decision_vector &, const double &);
		static void d_objfun_wrapper(const gsl_vector *, void *, gsl_vector *);
		static void fd_objfun_wrapper(const gsl_vector *, void *, double *, gsl_vector *);
		static void cleanup(gsl_vector *, gsl_multimin_fdfminimizer *);
		static void check_allocs(gsl_vector *, gsl_multimin_fdfminimizer *);
	private:
		const std::size_t	m_max_iter;
		const double		m_grad_tol;
		const double		m_numdiff_step_size;
		const double		m_step_size;
		const double		m_tol;
};

}}

#endif
