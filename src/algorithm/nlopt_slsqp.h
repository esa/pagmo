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

#ifndef PAGMO_ALGORITHM_NLOPT_SLSQP_H
#define PAGMO_ALGORITHM_NLOPT_SLSQP_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "base_nlopt.h"

namespace pagmo { namespace algorithm {

/// Wrapper for NLopt's Slsqp algorithm.
/**
 * From NLopt's documentation:
 *
 * <em>The algorithm optimizes successive second-order (quadratic/least-squares) approximations
 * of the objective function (via BFGS updates), with first-order (affine) approximations of the
 * constraints. The Fortran code was obtained from the SciPy project, who are responsible for
 * obtaining permission to distribute it under a free-software (3-clause BSD) license.
 * The code was modified for inclusion in NLopt by S. G. Johnson in 2010,  .... since roundoff
 * errors sometimes pushed SLSQP's parameters slightly outside the
 * bound constraints (not allowed by NLopt), we added checks to force the parameters within the
 * bounds. We fixed a bug in the LSEI subroutine (use of uninitialized variables) for the case
 * where the number of equality constraints equals the dimension of the problem. The LSQ
 * subroutine was modified to handle infinite lower/upper bounds (in which case those constraints
 * are omitted).</em>
 *
 * The inclusion in PaGMO required to write central difference code for the automated, numercal evaluation
 * of gradients. the rest was left unchaged
 *
 * NOTE: <em>Because the SLSQP code uses dense-matrix methods (ordinary BFGS, not low-storage BFGS),
 * it requires O(n2) storage and O(n3) time in n dimensions, which makes it less practical for 
 * optimizing more than a few thousand parameters.</em>
 *
 * This algorithm is a single-objective continuous minimiser that supports box constraints.
 *
 * @see http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#SLSQP
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */
class __PAGMO_VISIBLE nlopt_slsqp: public base_nlopt
{
	public:
		nlopt_slsqp(int = 100, const double & = 1E-6, const double & = 1E-6);
		base_ptr clone() const;
		std::string get_name() const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_nlopt>(*this);	
		}  
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::nlopt_slsqp)

#endif
