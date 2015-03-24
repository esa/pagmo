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

#ifndef PAGMO_ALGORITHM_GSL_BFGS2_H
#define PAGMO_ALGORITHM_GSL_BFGS2_H

#include <gsl/gsl_multimin.h>
#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "gsl_gradient.h"

namespace pagmo { namespace algorithm {

/// Wrapper for the GSL BFGS2 algorithm.
/**
 * @see algorithm::gsl_gradient for more information.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE gsl_bfgs2: public gsl_gradient
{
	public:
		gsl_bfgs2(int = 100, const double & = 1E-8, const double & = 1E-8, const double & = 0.01, const double & = 0.1);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		const gsl_multimin_fdfminimizer_type *get_gsl_minimiser_ptr() const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<gsl_gradient>(*this);
		}
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::gsl_bfgs2)

#endif
