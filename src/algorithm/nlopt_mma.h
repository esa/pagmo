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

#ifndef PAGMO_ALGORITHM_NLOPT_MMA_H
#define PAGMO_ALGORITHM_NLOPT_MMA_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "base_nlopt.h"

namespace pagmo { namespace algorithm {

/// Wrapper for NLopt's Method of Moving Asymptotes algorithm.
/**
 * From NLopt's documentation:
 *
 * <EM>This is an improved variant of the original MMA algorithm published by Svanberg in 1987,
 * which has become popular for topology optimization. (Note: "globally convergent" does not mean
 * that this algorithm converges to the global optimum; it means that it is guaranteed to converge
 * to some local minimum from any feasible starting point.) At each point x, MMA forms a local
 * approximation using the gradient of f and the constraint functions, plus a quadratic "penalty"
 * term to make the approximations "conservative" (upper bounds for the exact functions).
 * The precise approximation MMA forms is difficult to describe in a few words, because it
 * includes nonlinear terms consisting of a poles at some distance from x (outside of the current
 * trust region), almost a kind of Pade approximant. .</EM>
 *
 * The inclusion in PaGMO required to write central difference code for the automated, numercal evaluation
 * of gradients. the rest was left unchaged
 *
 * This algorithm is a single-objective continuous minimiser that supports box constraints.
 *
 * @see http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#MMA_.28Method_of_Moving_Asymptotes.29
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */
class __PAGMO_VISIBLE nlopt_mma: public base_nlopt
{
	public:
		nlopt_mma(int = 100, const double & = 1E-6, const double & = 1E-6);
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

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::nlopt_mma)

#endif
