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

#include <string>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <numeric>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <keplerian_toolbox/core_functions/propagate_lagrangian.h>
#include <keplerian_toolbox/astro_constants.h>
#include <keplerian_toolbox/lambert_problem.h>
#include <keplerian_toolbox/core_functions/array3D_operations.h>

#include "mga_target_event.h"


namespace pagmo { namespace problem {

mga_target_event::mga_target_event(
				 const kep_toolbox::plantes::planet_ptr start,
				 const kep_toolbox::planet::planet_ptr end,
				 const kep_toolbox::epoch t_end,
				 double T_max,
				 bool discount_launcher
		):base(6), m_start(start), m_end(end), m_t_end(t_end), m_T_max(T_max), m_discount_launcher(discount_launcher)
{
	// We check that all planets have equal central body
	if (start->get_mu_central_body() != end->get_mu_central_body()) {
		pagmo_throw(value_error,"The planets do not all have the same mu_central_body");  
	}

	// We check that T_max is positive
	if (T_max <= 0) {
		pagmo_throw(value_error,"T_max must be larger than zero");
	}

	// Now setting the problem bounds
	decision_vector lb(6,0.0), ub(6,1.0);
	lb[0] = 1.; lb[5] = 1e-5;
	ub[0] = T_max; ub[2] = 12000.0; ub[5] = 1-1e-5;

	set_bounds(lb,ub);
}

/// Copy Constructor. Performs a deep copy
mga_target_event::mga_target_event(const mga_target_event &p) : base(p.get_dimension()), m_t_end(p.m_t_end), m_T_max(p.m_T_max), m_discount_launcher(p.m_discount_launcher)
{
	set_bounds(p.get_lb(),p.get_ub());
	m_start = p.m_start->clone();
	m_end = p.m_end->clone();
}

/// Clone method.
base_ptr mga_target_event::clone() const
{
	return base_ptr(new mga_target_event(*this));
}

/// Implementation of the objective function.
void mga_target_event::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	//1 -  we 'decode' the chromosome
	double tw = x[0]*x[1];				//transforms the waiting time from a fraction to days

	//2 - We compute the epochs and ephemerides of the planetary encounters
	kep_toolbox::epoch t0 = kep_toolbox::epoch(m_t_end.mjd2000() - x[0]);	//epoch of departure
	double tof = (x[0] - tw);												//in days
//std::cout << "tof: " << tof << std::endl;
	kep_toolbox::epoch t1 = kep_toolbox::epoch(t0.mjd2000() + tof);			//epoch

	kep_toolbox::array3D r0,v0,r1,v1,v_inf_vett,v0_sc;
	m_start->eph(t0,r0,v0);
	m_end->eph(t1,r1,v1);
//std::cout << "r0: " << r0 << std::endl;
//std::cout << "v0: " << v0 << std::endl;
//std::cout << "r1: " << r1 << std::endl;
//std::cout << "v1: " << v1 << std::endl;

	//3 - We compute the absolute velocity at departure
	double theta = 2*M_PI*x[3];
	double phi = acos(2. * x[4]-1.) - M_PI / 2.;
	v_inf_vett[0] = x[2]*cos(phi)*cos(theta);
	v_inf_vett[1] = x[2]*cos(phi)*sin(theta);
	v_inf_vett[2] = x[2]*sin(phi);
	kep_toolbox::sum(v0_sc,v0,v_inf_vett);

//std::cout << "v_inf_vett: " << v_inf_vett << std::endl;
//std::cout << "v_sc: " << v0_sc << std::endl;
//std::cout << "tof: " << tof << std::endl;
//std::cout << "x[5]: " << x[5] << std::endl;
//std::cout << "x: " << x << std::endl;

	//4 - And we propagate up to the DSM positon to then solve a Lambert problem
	kep_toolbox::propagate_lagrangian(r0,v0_sc,tof*x[5]*ASTRO_DAY2SEC, ASTRO_MU_SUN);
	kep_toolbox::lambert_problem l(r0,r1,tof*(1-x[5])*ASTRO_DAY2SEC,ASTRO_MU_SUN, false,0);

//std::cout << "r0: " << r0 << std::endl;
//std::cout << "v0: " << v0_sc << std::endl;

	double DV1 = x[2];
	kep_toolbox::array3D dv2,dv3;
	kep_toolbox::diff(dv2,v0_sc,l.get_v1()[0]);
	double DV2 = kep_toolbox::norm(dv2);
	kep_toolbox::diff (dv3,v1,l.get_v2()[0]);
	double DV3 = kep_toolbox::norm(dv3);
	if (m_discount_launcher) {
		DV1 = std::max(0.,DV1-6000.);
	}
	f[0] = DV1+DV2+DV3;
//	std::cout << "DV1: " << DV1 << std::endl;
//	std::cout << "DV2: " << DV2 << std::endl;
//	std::cout << "DV3: " << DV3 << std::endl;
}

std::string mga_target_event::get_name() const
{
	return "MGA-Target-Event";
}


/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the values vector, the weights vectors and the max weight. It is concatenated
 * with the base::problem human_readable
 */
std::string mga_target_event::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "\n\tStart: " << m_start->get_name();
	oss << "\n\tEnd: " + m_end->get_name();
	oss << "\n\tFinal epoch: " << boost::lexical_cast<std::string>(m_t_end.mjd2000());
	oss << "\n\tDiscount launcher dv?: " << (m_discount_launcher ? "True" : "False");
	return oss.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::mga_target_event)

