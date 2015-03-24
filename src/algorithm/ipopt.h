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

#ifndef PAGMO_ALGORITHM_IPOPT_H
#define PAGMO_ALGORITHM_IPOPT_H

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"

#include <coin/IpIpoptApplication.hpp>



namespace pagmo { namespace algorithm {

/// Wrapper for the IPOPT solver
/**
 * IPOPT is an Interior Point Optimization solver released under the CPL licence. Our wrapper is
 * generic and does not use the hessian information letting IPOPT deal with it. Gradients and Jacobian are
 * evaluated by central differences with a stepsize of 10-8 stepsize. We offer control
 * over only a few IPOPT options, take care in particular of the dual_inf_tol as the hessian_approximation,
 * limited_memory option (i.e. hessian is not provided by the user) makes the dual_inf to remain quite high.
 *
 * CAREFUL:
 *	1 - IPOPT works only for minimization.
 *	2 - The final solution is guaranteed to be within the box constraints by forcing it after the ipopt call
 *
 * IPOPT has severa lib dependencies that overlaps with SNOPT lib dependencies. Take care that you link
 * to the correct ones, especially for the blas lib that can give problem if coming from SNOPT native source code
 *
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 *
 */

class __PAGMO_VISIBLE ipopt: public base
{
public:

	ipopt(const int &max_iter = 100, 
		const double &constr_viol_tol = 1e-8, 
		const double &dual_inf_tol = 1e-8, 
		const double &compl_inf_tol = 1e-8,
		const bool &nlp_scaling_method = true,
		const double &obj_scaling_factor = 1.0,
		const double &mu_init = 0.1
 	    );
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;

protected:
	std::string human_readable_extra() const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_max_iter);
		ar & const_cast<double &>(m_constr_viol_tol);
		ar & const_cast<double &>(m_dual_inf_tol);
		ar & const_cast<double &>(m_compl_inf_tol);
		ar & const_cast<bool &>(m_nlp_scaling_method);
		ar & const_cast<double &>(m_obj_scaling_factor);
		ar & const_cast<double &>(m_mu_init);
	}  
	const int m_max_iter;
	const double m_constr_viol_tol;
	const double m_dual_inf_tol;
	const double m_compl_inf_tol;
	const bool m_nlp_scaling_method;
	const double m_obj_scaling_factor;
	const double m_mu_init;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::ipopt)

#endif
