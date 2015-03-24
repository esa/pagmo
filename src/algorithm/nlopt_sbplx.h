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

#ifndef PAGMO_ALGORITHM_NLOPT_SBPLX_H
#define PAGMO_ALGORITHM_NLOPT_SBPLX_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "base_nlopt.h"

namespace pagmo { namespace algorithm {

/// Wrapper for NLopt's Sbplx algorithm.
/**
 * From NLopt's documentation:
 *
 * <EM>This is my re-implementation of Tom Rowan's "Subplex" algorithm. As Rowan expressed a preference that other implementations of
 * his algorithm use a different name, I called my implementation "Sbplx" [...].
 * Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces) is claimed to be much more efficient and robust
 * than the original Nelder-Mead, while retaining the latter's facility with discontinuous objectives, and in my experience these claims
 * seem to be true in many cases.</EM>
 *
 * This algorithm is a derivative-free single-objective continuous minimiser that supports box constraints.
 *
 * @see T. Rowan, "Functional Stability Analysis of Numerical Algorithms", Ph.D. thesis, Department of Computer Sciences, University of Texas at Austin, 1990.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE nlopt_sbplx: public base_nlopt
{
	public:
		nlopt_sbplx(int = 100, const double & = 1E-6, const double & = 1E-6);
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

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::nlopt_sbplx)

#endif
