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

#ifndef PAGMO_ALGORITHM_NLOPT_AUG_LAG_EQ_H
#define PAGMO_ALGORITHM_NLOPT_AUG_LAG_EQ_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "base_nlopt.h"

namespace pagmo { namespace algorithm {

/// Wrapper for NLopt's Augmented Lagrangian algorithm (using penalties only for the equalities)
/**
 * From NLopt's documentation:
 *
 * <EM>This method combines the objective function and the nonlinear inequality/equality constraints (if any)
 * in to a single function: essentially, the objective plus a "penalty" for any violated constraints.
 * This modified objective function is then passed to another optimization algorithm with no nonlinear
 * constraints [the auxiliary algorithm]. If the constraints are violated by the solution of this sub-problem, then the size
 * of the penalties is increased and the process is repeated; eventually, the process must converge
 * to the desired solution (if it exists).</EM>
 *
 * The inclusion in PaGMO required to give a fixed choice forthe auxiliary algorithms as nlopt interface is
 * too different from pagmo's to allow pagmo algorithm being passed as auxiliaries.
 *
 * This version of the augmented lagrangian delegates to the auxiliary algorithm the inequality contraints handling
 *
 * This algorithm is a single-objective continuous minimiser that supports any type of constraints
 *
 * @see http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Augmented_Lagrangian_algorithm
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */
class __PAGMO_VISIBLE nlopt_aug_lag_eq: public base_nlopt
{
	public:
		nlopt_aug_lag_eq(int=1, int = 100, const double & = 1E-6, const double & = 1E-6, int = 100, const double & = 1E-6, const double & = 1E-6);
		base_ptr clone() const;
		std::string get_name() const;
		void set_local(size_t) const;
		std::string human_readable_extra() const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_nlopt>(*this);
			ar & m_aux_algo_id;
			ar & const_cast<std::size_t &>(m_aux_max_iter);
			ar & const_cast<double &>(m_aux_ftol);
			ar & const_cast<double &>(m_aux_xtol);
		}  
		int m_aux_algo_id;
		const std::size_t	m_aux_max_iter;
		const double		m_aux_ftol;
		const double		m_aux_xtol;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::nlopt_aug_lag_eq)

#endif
