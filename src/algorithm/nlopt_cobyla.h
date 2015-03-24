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

#ifndef PAGMO_ALGORITHM_NLOPT_COBYLA_H
#define PAGMO_ALGORITHM_NLOPT_COBYLA_H

#include "../config.h"
#include "../serialization.h"
#include "base_nlopt.h"

namespace pagmo { namespace algorithm {

/// Wrapper for NLopt's COBYLA algorithm.
/**
 * COBYLA (Constrained Optimization BY Linear Approximations) is a derivative-free single-objective continuous minimiser that supports nonlinear inequality constraints.
 *
 * @see M. J. D. Powell, "A direct search optimization method that models the objective and constraint functions by linear interpolation," in Advances in Optimization and Numerical Analysis, eds. S. Gomez and J.-P. Hennart (Kluwer Academic: Dordrecht, 1994), p. 51-67.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE nlopt_cobyla: public base_nlopt
{
	public:
		nlopt_cobyla(int = 100, const double & = 1E-6, const double & = 1E-6);
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

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::nlopt_cobyla)

#endif
