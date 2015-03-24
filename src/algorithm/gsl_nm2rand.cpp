/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#include <gsl/gsl_multimin.h>

#include "../exceptions.h"
#include "gsl_derivative_free.h"
#include "gsl_nm2rand.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Will invoke internally the constructor from algorithm::gsl_derivative_free with the specified parameters.
 *
 * @see gsl_gradient::gsl_derivative_free().
 */
gsl_nm2rand::gsl_nm2rand(int max_iter, const double &tol, const double &step_size):
	gsl_derivative_free(max_iter,tol,step_size) {}

/// Clone method.
/**
 * @return algorithm::base_ptr to a copy of this.
 */
base_ptr gsl_nm2rand::clone() const
{
	return base_ptr(new gsl_nm2rand(*this));
}

const gsl_multimin_fminimizer_type *gsl_nm2rand::get_gsl_minimiser_ptr() const
{
	return gsl_multimin_fminimizer_nmsimplex2rand;
}

/// Algorithm name
std::string gsl_nm2rand::get_name() const
{
	return "Nelder-Mead simplex (GSL, version 2 + randomly-oriented initial simplex)";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::gsl_nm2rand)
