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

#ifndef PAGMO_ALGORITHM_GSL_NM2_H
#define PAGMO_ALGORITHM_GSL_NM2_H

#include <gsl/gsl_multimin.h>

#include "../config.h"
#include "../serialization.h"
#include "gsl_derivative_free.h"

namespace pagmo { namespace algorithm {

/// Wrapper for the GSL Nelder-Mead simplex algorithm (version 2).
/**
 * @see algorithm::gsl_derivative_free for more information.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE gsl_nm2: public gsl_derivative_free
{
	public:
		gsl_nm2(int max_iter = 100, const double &tol = 1E-6, const double &step_size = 1);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		const gsl_multimin_fminimizer_type *get_gsl_minimiser_ptr() const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<gsl_derivative_free>(*this);
		}
};

}}
BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::gsl_nm2)

#endif
