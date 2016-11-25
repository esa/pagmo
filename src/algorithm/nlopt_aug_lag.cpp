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

#include <nlopt.hpp>

#include "base_nlopt.h"
#include "nlopt_aug_lag.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify some of the parameters of the IPOPT solver. The algorithm stopping criteria will
 * be 1) if the number of iteration exceeds max_iter 2) the three tolerances are met
 *
 * @param[in] aux_algo_id Auxiliary algorithm (1=SPBLX, 2=COBYLA, 3=BOBYQA, 4=LBFGS)
 * @param[in] max_iter stop-criteria (number of iterations)
 * @param[in] ftol stop-criteria (absolute on the obj-fun)
 * @param[in] xtol stop-criteria (absolute on the chromosome)
 * @param[in] aux_max_iter stop-criteria for the auxiliary algorithm (number of iterations)
 * @param[in] aux_ftol stop-criteria for the auxiliary algorithm (number of iterations)
 * @param[in] aux_xtol stop-criteria for the auxiliary algorithm (number of iterations)
 * @throws value_error if max_iter or tolerances are negative
 *
 * @see pagmo::algorithm::base_nlopt::base_nlopt()
 */
nlopt_aug_lag::nlopt_aug_lag(int aux_algo_id, int max_iter, const double &ftol, const double &xtol, int aux_max_iter, const double &aux_ftol, const double &aux_xtol):base_nlopt(nlopt::AUGLAG,true,false,max_iter,ftol,xtol), m_aux_algo_id(aux_algo_id), m_aux_max_iter(aux_max_iter), m_aux_ftol(aux_ftol), m_aux_xtol(aux_xtol) {
	if ( (aux_ftol <= 0) || (aux_xtol <= 0) ) {
		pagmo_throw(value_error,"tolerances for the local optimizer must be positive");
	}
	if ((aux_algo_id >4)||(aux_algo_id<1)) {
		pagmo_throw(value_error,"local algorithm id must be one of 1.2.3.4");
	}
}

base_ptr nlopt_aug_lag::clone() const
{
	return base_ptr(new nlopt_aug_lag(*this));
}

/// Set the local optimizer
void nlopt_aug_lag::set_local(size_t d) const{
	nlopt::opt aux_opt(nlopt::LN_SBPLX,1);
	switch(m_aux_algo_id)
	{
		case 1:
			aux_opt = nlopt::opt(nlopt::LN_SBPLX,d);
		break;
		case 2:
			aux_opt= nlopt::opt(nlopt::LN_COBYLA,d);
		break;
		case 3:
			aux_opt = nlopt::opt(nlopt::LN_BOBYQA,d);
		break;
		case 4:
			aux_opt = nlopt::opt(nlopt::LD_LBFGS,d);
		break;
	}
	aux_opt.set_ftol_abs(m_aux_ftol);
	aux_opt.set_xtol_abs(m_aux_xtol);
	aux_opt.set_maxeval(m_aux_max_iter);
	m_opt.set_local_optimizer(aux_opt);
}

std::string nlopt_aug_lag::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "max_iter:" << m_max_iter << ' ';
	oss << "ftol:" << m_ftol << " ";
	oss << "xtol:" << m_xtol << " ";
	oss << "aux_max_iter:" << m_aux_max_iter << ' ';
	oss << "aux_ftol:" << m_aux_ftol << " ";
	oss << "aux_xtol:" << m_aux_xtol;

	return oss.str();
}

/// Algorithm name
std::string nlopt_aug_lag::get_name() const
{
	std::string local;
	switch(m_aux_algo_id)
	{
		case 1:
			local = "sbplx";
		break;
		case 2:
			local = "cobyla";
		break;
		case 3:
			local = "bobyqa";
		break;
		case 4:
			local = "lbfgs";
		break;
	}
	
	return "Augmented Lagrangian - " + local + " (NLOPT)";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::nlopt_aug_lag)
