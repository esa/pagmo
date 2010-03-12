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

#ifndef PAGMO_ALGORITHM_BASE_GSL_H
#define PAGMO_ALGORITHM_BASE_GSL_H

#include <gsl/gsl_vector.h>

#include "../gsl_init.h"
#include "../problem/base.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Base class for GSL algorithms.
/**
 * All GSL algorithms should derive from this class, which will automatically setup
 * the GSL environment for use in PaGMO, and will provide building blocks for wrapping the GSL minimisers.
 * Please note that GSL provides function minimisers, so that regardless of the comparison functions implemented in
 * problem::base, all GSL algorithms will try to minimise the objective function.
 *
 * This class of algorithms supports single-objective, unconstrained, continuous optimisation.
 *
 * @see http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html for an overview of the minimisers available in the GSL.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class base_gsl: public base
{
	public:
		base_gsl();
	protected:
		/// Structure to feed parameters to the wrappers for the objective function and its derivative.
		struct objfun_wrapper_params
		{
			/// Pointer to the problem.
			problem::base const		*p;
			/// Decision vector.
			decision_vector			x;
			/// Fitness vector.
			fitness_vector			f;
			/// Step size for the computation of the gradient.
			double				step_size;
		};
		static double objfun_wrapper(const gsl_vector *, void *);
	private:
		static gsl_init g_init;
};

}}

#endif
