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

#include "base_stochastic.h"
#include "../serialization.h"

namespace pagmo { namespace problem {


/// Constructor from global dimension and random seed
/**
 * Lower and upper bounds are set to 0 and 1 respectively. 
 * The problem built is unconstrained, single objective and with no integer dimension. 
 *
 * @param[in] dim global dimension of the problem.
 * @param[in] seed random number generator seed
 */
base_stochastic::base_stochastic(int dim, unsigned int seed) : base(dim), m_drng(seed), m_urng(seed), m_seed(seed)
{
}

/// Constructor from global dimension, integer dimension, fitness dimension, global constraints dimension, inequality constraints dimension and constraints tolerance.
/**
 * @param[in] n global dimension of the problem.
 * @param[in] ni dimension of the combinatorial part of the problem.
 * @param[in] nf dimension of the fitness vector of the problem.
 * @param[in] nc global number of constraints.
 * @param[in] nic number of inequality constraints.
 * @param[in] c_tol constraints tolerance (equal for all constraints)
 * @param[in] seed random number generator seed
*/
base_stochastic::base_stochastic(int n, int ni, int nf, int nc, int nic, const double &c_tol, unsigned int seed): base((int)n, ni, nf, nc, nic, c_tol), m_drng(seed), m_urng(seed), m_seed(seed)
{
}

/// Constructor from global dimension, integer dimension, fitness dimension, global constraints dimension, inequality constraints dimension and constraints tolerance vector
/**
 * @param[in] n global dimension of the problem.
 * @param[in] ni dimension of the combinatorial part of the problem.
 * @param[in] nf dimension of the fitness vector of the problem.
 * @param[in] nc global number of constraints.
 * @param[in] nic number of inequality constraints.
 * @param[in] c_tol constraints tolerance std::vector
 * @param[in] seed random number generator seed
*/
base_stochastic::base_stochastic(int n, int ni, int nf, int nc, int nic, const std::vector<double> &c_tol, unsigned int seed): base((int)n, ni, nf, nc, nic, c_tol), m_drng(seed), m_urng(seed), m_seed(seed)
{
}

/// Sets the pseudo random generator seed
/**
 * Sets the pseudo random generator seed.
 *
 * @param[in] seed random number generator seed
 */
void base_stochastic::set_seed(unsigned int seed) const {
	m_seed = seed;
	// As the problem is now muted we must reset the caches that contain the evaluations w.r.t. the old seed
	reset_caches();
}

/// Gets the pseudo random generator seed
/**
 * Gets the pseudo random generator seed.
 *
 * @return random number generator seed
 */
unsigned int base_stochastic::get_seed() const {
	return m_seed;
}

}} //namespaces

//BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::base_stochastic);
