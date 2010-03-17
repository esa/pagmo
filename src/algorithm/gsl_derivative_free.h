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

#ifndef PAGMO_ALGORITHM_GSL_DERIVATIVE_FREE_H
#define PAGMO_ALGORITHM_GSL_DERIVATIVE_FREE_H

#include <cstddef>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <string>

#include "../config.h"
#include "../population.h"
#include "base_gsl.h"

namespace pagmo { namespace algorithm {

/// Wrapper for GSL minimisers without derivatives.
/**
 * This class can be used to build easily a wrapper around a GSL minimiser without derivatives.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE gsl_derivative_free: public base_gsl
{
	protected:
		gsl_derivative_free(const gsl_multimin_fminimizer_type *, int, const double &, const double &);
		std::string human_readable_extra() const;
		void evolve(population &) const;
	private:
		static void cleanup(gsl_vector *, gsl_vector *, gsl_multimin_fminimizer *);
		static void check_allocs(gsl_vector *, gsl_vector *, gsl_multimin_fminimizer *);
	private:
		const gsl_multimin_fminimizer_type	*m_minimiser;
		const std::size_t			m_max_iter;
		const double				m_tol;
		const double				m_step_size;
};

} }

#endif
