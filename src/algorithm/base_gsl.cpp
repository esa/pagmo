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

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <gsl/gsl_vector.h>
#include <stdexcept>

#include "../exceptions.h"
#include "../gsl_init.h"
#include "../problem/base.h"
#include "base_gsl.h"
#include "base.h"

namespace pagmo { namespace algorithm {

// Static mutex.
boost::mutex base_gsl::m_mutex;

const gsl_init &base_gsl::get_init()
{
	// Protect from simultaneous access.
	boost::lock_guard<boost::mutex> lock(m_mutex);
	static const gsl_init retval;
	return retval;
}

/// Default constructor.
/**
 * Will make sure that the GSL environment is initialised properly.
 */
base_gsl::base_gsl():base()
{
	if (!get_init().m_init) {
		pagmo_throw(std::runtime_error,"GSL support could not be initialised");
	}
}

/// Objective function wrapper.
/**
 * @param[in] v pointer to the gsl_vector representing the decision vector.
 * @param[in] params pointer to extra parameters for the internal function.
 *
 * @return the fitness of the input decision vector.
 */
double base_gsl::objfun_wrapper(const gsl_vector *v, void *params)
{
	objfun_wrapper_params *par = (objfun_wrapper_params *)params;
	// Size of the continuous part of the problem.
	const problem::base::size_type cont_size = par->p->get_dimension() - par->p->get_i_dimension();
	// Fill up the continuous part of temporary storage with the contents of v.
	for (problem::base::size_type i = 0; i < cont_size; ++i) {
		par->x[i] = gsl_vector_get(v,i);
	}
	// Compute the objective function.
	par->p->objfun(par->f,par->x);
	return par->f[0];
}

}}
