/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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
 * @see pagmo::algorithm::base_nlopt::base_nlopt()
 */
nlopt_aug_lag::nlopt_aug_lag(int l_algo_id, int max_iter, const double &ftol, const double &xtol, int l_max_iter, const double &l_ftol, const double &l_xtol):base_nlopt(nlopt::AUGLAG,true,false,max_iter,ftol,xtol), m_l_algo_id(l_algo_id), m_l_max_iter(l_max_iter), m_l_ftol(l_ftol), m_l_xtol(l_xtol) {
	if ( (l_ftol <= 0) || (l_xtol <= 0) ) {
		pagmo_throw(value_error,"tolerances for the local optimizer must be positive");
	}
	if ((l_algo_id >4)||(l_algo_id<1)) {
		pagmo_throw(value_error,"local algorithm id must be one of 1.2.3.4");
	}
}

base_ptr nlopt_aug_lag::clone() const
{
	return base_ptr(new nlopt_aug_lag(*this));
}

/// Set the local optimizer
void nlopt_aug_lag::set_local(size_t d) const{
	nlopt::opt l_opt(nlopt::LN_SBPLX,1);
	switch(m_l_algo_id)
	{
		case 1:
			l_opt = nlopt::opt(nlopt::LN_SBPLX,d);
		break;
		case 2:
			l_opt= nlopt::opt(nlopt::LN_COBYLA,d);
		break;
		case 3:
			l_opt = nlopt::opt(nlopt::LN_BOBYQA,d);
		break;
		case 4:
			l_opt = nlopt::opt(nlopt::LD_LBFGS,d);
		break;
	}
	l_opt.set_ftol_abs(m_l_ftol);
	l_opt.set_xtol_abs(m_l_xtol);
	l_opt.set_maxeval(m_l_max_iter);
	m_opt.set_local_optimizer(l_opt);
}

/// Algorithm name
std::string nlopt_aug_lag::get_name() const
{
	std::string local;
	switch(m_l_algo_id)
	{
		case 1:
			local = "sbplx";
		break;
		case 2:
			local = "cobyla";
		break;
		case 3:
			local = "bobya";
		break;
		case 4:
			local = "lbfgs";
		break;
	}
	
	return "Augented Lagrangian - " + local + " (NLOPT)";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::nlopt_aug_lag);
