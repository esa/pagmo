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
#include <keplerian_toolbox/core_functions/propagate_lagrangian.h>
#include <keplerian_toolbox/core_functions/fb_prop.h>
#include <keplerian_toolbox/core_functions/closest_distance.h>
#include <keplerian_toolbox/lambert_problem.h>

#include "mga_incipit.h"


namespace pagmo { namespace problem {

 
/// Constructor
/**
 * Constructs a global optimization problem (box-bounded, continuous) representing an interplanetary trajectory modelled
 * as a Multiple Gravity Assist trajectory that allows one only Deep Space Manouvre between each leg.
 *  
 * @param[in] seq std::vector of kep_toolbox::planet_ptr containing the encounter sequence for the trajectoty (including the initial planet)
 * @param[in] t0_l kep_toolbox::epoch representing the lower bound for the launch epoch
 * @param[in] t0_u kep_toolbox::epoch representing the upper bound for the launch epoch
 * @param[in] tof time-of-flight vector containing lower and upper bounds (in days) for the various legs time of flights
 * 
 * @throws value_error if the planets in seq do not all have the same central body gravitational constant
 * @throws value_error if tof has a size different from seq.size()
 */
mga_incipit::mga_incipit(const std::vector<kep_toolbox::planet::planet_ptr> seq, 
			 const kep_toolbox::epoch t0_l, const kep_toolbox::epoch t0_u,
			 const std::vector<std::vector<double> > tof
			): 
			 base(4*seq.size()), m_tof(tof)
{
	// We check that all planets have equal central body
	std::vector<double> mus(seq.size());
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i< seq.size(); ++i) {
		mus[i] = seq[i]->get_mu_central_body();
	}
	if ((unsigned int)std::count(mus.begin(), mus.end(), mus[0]) != mus.size()) {
		pagmo_throw(value_error,"The planets do not all have the same mu_central_body");  
	}
	// We check the consistency of the time of flights
	if (tof.size() != seq.size()) {
		pagmo_throw(value_error,"The time-of-flight vector (tof) has the wrong length");  
	}
	for (size_t i = 0; i < tof.size(); ++i) {
		if (tof[i].size()!=2) pagmo_throw(value_error,"Each element of the time-of-flight vector (tof)  needs to have dimension 2 (lower and upper bound)"); 
	}
	
	// Filling in the planetary sequence data member. This requires to construct the polymorphic planets via their clone method 
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < seq.size(); ++i) {
		m_seq.push_back(seq[i]->clone());
	}
	
	// Now setting the problem bounds
	size_type dim(4*m_tof.size());
	decision_vector lb(dim), ub(dim);
	
	// First leg
	lb[0] = t0_l.mjd2000(); ub[0] = t0_u.mjd2000();
	lb[1] = 0; lb[2] = 0; ub[1] = 1; ub[2] = 1;
	lb[3] = m_tof[0][0]; ub[3] = m_tof[0][1];
	
	// Successive legs
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 1; i < m_tof.size(); ++i) {
		lb[4*i] = - 2 * boost::math::constants::pi<double>();    ub[4*i] = 2 * boost::math::constants::pi<double>();
		lb[4*i+1] = 1.1;  ub[4*i+1] = 30;
		lb[4*i+2] = 1e-5; ub[4*i+2] = 1-1e-5;
		lb[4*i+3] = m_tof[i][0]; ub[4*i+3] = m_tof[i][1];
	}
	
	// Adjusting the minimum and maximum allowed fly-by rp to the one defined in the kep_toolbox::planet class
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < m_tof.size()-1; ++i) {
		lb[4*i+5] = m_seq[i]->get_safe_radius() / m_seq[i]->get_radius();
		ub[4*i+5] = (m_seq[i]->get_radius() + 2000000) / m_seq[i]->get_radius(); //from gtoc6 problem description
	}
	set_bounds(lb,ub);
}

/// Copy Constructor. Performs a deep copy
mga_incipit::mga_incipit(const mga_incipit &p) : base(p.get_dimension()), m_tof(p.m_tof)
{
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < p.m_seq.size();++i) {
		m_seq.push_back(p.m_seq[i]->clone());
	}
	set_bounds(p.get_lb(),p.get_ub());
}

/// Clone method.
base_ptr mga_incipit::clone() const
{
	return base_ptr(new mga_incipit(*this));
}

/// Implementation of the objective function.
void mga_incipit::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
try {
	double common_mu = m_seq[0]->get_mu_central_body();
	// 1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
	std::vector<double> T(m_seq.size(),0.0);
	
	for (size_t i = 0; i<m_seq.size(); ++i) {
		T[i] = x[4*i+3];
	}
	// 2 - We compute the epochs and ephemerides of the planetary encounters
	std::vector<kep_toolbox::epoch>   t_P(m_seq.size());
	std::vector<kep_toolbox::array3D> r_P(m_seq.size());
	std::vector<kep_toolbox::array3D> v_P(m_seq.size());
	std::vector<double> DV(m_seq.size());
	for (size_t i = 0; i<r_P.size(); ++i) {
		t_P[i] = kep_toolbox::epoch(x[0] + std::accumulate(T.begin(),T.begin()+1+i,0.0));
		m_seq[i]->eph(t_P[i], r_P[i], v_P[i]);
	}

	// 3 - We start with the first leg
	double theta = 2*boost::math::constants::pi<double>()*x[1];
	double phi = acos(2*x[2]-1)-boost::math::constants::pi<double>() / 2;
	double d,d2,ra,ra2;
	kep_toolbox::array3D r = { {ASTRO_JR*1000*cos(phi)*sin(theta), ASTRO_JR*1000*cos(phi)*cos(theta), ASTRO_JR*1000*sin(phi)} };
	kep_toolbox::array3D v;
	kep_toolbox::lambert_problem l(r,r_P[0],T[0]*ASTRO_DAY2SEC,common_mu,false,false);
	kep_toolbox::array3D v_beg_l = l.get_v1()[0];
	kep_toolbox::array3D v_end_l = l.get_v2()[0];

	DV[0] = std::abs(kep_toolbox::norm(v_beg_l)-3400.0);
	
	// 4 - And we proceed with each successive leg (if any)
	kep_toolbox::array3D v_out;
	for (size_t i = 1; i<m_seq.size(); ++i) {
		// Fly-by
		kep_toolbox::fb_prop(v_out, v_end_l, v_P[i-1], x[4*i+1] * m_seq[i-1]->get_radius(), x[4*i], m_seq[i-1]->get_mu_self());
	    r = r_P[i-1];
		v = v_out;
		// s/c propagation before the DSM
		kep_toolbox::propagate_lagrangian(r,v,x[4*i+2]*T[i]*ASTRO_DAY2SEC,common_mu);
		kep_toolbox::closest_distance(d, ra, r_P[i-1], v_out, r, v, common_mu);

		// Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
		double dt = (1-x[4*i+2])*T[i]*ASTRO_DAY2SEC;
		kep_toolbox::lambert_problem l2(r,r_P[i],dt,common_mu,false,false);
		v_end_l = l2.get_v2()[0];
		v_beg_l = l2.get_v1()[0];
		kep_toolbox::closest_distance(d2,ra2,r,v_beg_l, r_P[i], v_end_l, common_mu);
		if (d < d2)
		{
			d = d/ASTRO_JR;
		} else {
			d = d2/ASTRO_JR;
		}

		// DSM occuring at time nu2*T2
		kep_toolbox::diff(v_out, v_beg_l, v);
		DV[i] = kep_toolbox::norm(v_out) + std::max((2.0-d),0.0) * 1000.0;
	}
	// Now we return the objective(s) function
	f[0] = std::accumulate(DV.begin(),DV.end(),0.0); 
//Here the lambert solver or the lagrangian propagator went wrong
} catch (...) {
	f[0] = boost::numeric::bounds<double>::highest();
} 
}

/// Outputs a stream with the trajectory data
/**
 * While the chromosome contains all necessary information to describe a trajectory, mission analysits
 * often require a different set of data to evaluate its use. This method outputs a stream with
 * information on the trajectory that is otherwise 'hidden' in the chromosome
 *
 * \param[in] x chromosome representing the trajectory in the optimization process
 * \returns an std::string with launch dates, DV magnitues and other information on the trajectory
 */

std::string mga_incipit::pretty(const std::vector<double> &x) const {
  
	// We set the std output format
	std::ostringstream s;
	s.precision(15);
	s << std::scientific;
	
	double d,ra,d2,ra2;

	double common_mu = m_seq[0]->get_mu_central_body();
	// 1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
	std::vector<double> T(m_seq.size(),0.0);
	
	for (size_t i = 0; i<m_seq.size(); ++i) {
		T[i] = x[4*i+3];
	}
	// 2 - We compute the epochs and ephemerides of the planetary encounters
	std::vector<kep_toolbox::epoch>   t_P(m_seq.size());
	std::vector<kep_toolbox::array3D> r_P(m_seq.size());
	std::vector<kep_toolbox::array3D> v_P(m_seq.size());
	std::vector<double> DV(m_seq.size());
	for (size_t i = 0; i<r_P.size(); ++i) {
		t_P[i] = kep_toolbox::epoch(x[0] + std::accumulate(T.begin(),T.begin()+1+i,0.0));
		m_seq[i]->eph(t_P[i], r_P[i], v_P[i]);
	}

	// 3 - We start with the first leg
	double theta = 2*boost::math::constants::pi<double>()*x[1];
	double phi = acos(2*x[2]-1)-boost::math::constants::pi<double>() / 2;
	kep_toolbox::array3D r = { {ASTRO_JR * 1000*cos(phi)*sin(theta), ASTRO_JR * 1000*cos(phi)*cos(theta), ASTRO_JR * 1000*sin(phi)} };
	kep_toolbox::array3D v;
	
	kep_toolbox::lambert_problem l(r,r_P[0],T[0]*ASTRO_DAY2SEC,common_mu,false,false);
	kep_toolbox::array3D v_beg_l = l.get_v1()[0];
	kep_toolbox::array3D v_end_l = l.get_v2()[0];
	kep_toolbox::closest_distance(d,ra,r,v_beg_l, r_P[0], v_end_l, common_mu);

	DV[0] = std::abs(kep_toolbox::norm(v_beg_l)-3400.0);
	kep_toolbox::array3D v_out,mem_vin,mem_vout,mem_vP;
	
	s << "\nFirst Leg: 1000JR to " << m_seq[0]->get_name() << std::endl; 
	s << "\tDeparture: " << kep_toolbox::epoch(x[0]) << " (" << x[0] << " mjd2000) " << std::endl; 
	s << "\tDuration: " << T[0] << "days" << std::endl; 
	s << "\tInitial Velocity Increment (m/s): " << DV[0] << std::endl; 
	kep_toolbox::diff(v_out, v_end_l, v_P[0]);
	s << "\tArrival relative velocity at " << m_seq[0]->get_name() << " (m/s): " << kep_toolbox::norm(v_out)  << std::endl; 
	s << "\tClosest distance: " << d / ASTRO_JR;
	
	// 4 - And we proceed with each successive leg (if any)
	for (size_t i = 1; i<m_seq.size(); ++i) {
		// Fly-by
		kep_toolbox::fb_prop(v_out, v_end_l, v_P[i-1], x[4*i+1] * m_seq[i-1]->get_radius(), x[4*i], m_seq[i-1]->get_mu_self());
		// s/c propagation before the DSM
		r = r_P[i-1];
		v = v_out;
		mem_vout = v_out;
		mem_vin = v_end_l;
		mem_vP = v_P[i-1];
		
		kep_toolbox::propagate_lagrangian(r,v,x[4*i+2]*T[i]*ASTRO_DAY2SEC,common_mu);
		kep_toolbox::closest_distance(d, ra, r_P[i-1], v_out, r, v, common_mu);
		
		// Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
		double dt = (1-x[4*i+2])*T[i]*ASTRO_DAY2SEC;
		kep_toolbox::lambert_problem l2(r,r_P[i],dt,common_mu,false,false);
		v_end_l = l2.get_v2()[0];
		v_beg_l = l2.get_v1()[0];
		kep_toolbox::closest_distance(d2,ra2,r,v_beg_l, r_P[i], v_end_l, common_mu);
		
		if (d < d2)
		{
			d = d/ASTRO_JR;
			ra = ra/ASTRO_JR;
		} else {
			d = d2/ASTRO_JR;
			ra = ra2/ASTRO_JR;
		}
		// DSM occuring at time nu2*T2
		kep_toolbox::diff(v_out, v_beg_l, v);
		DV[i] = kep_toolbox::norm(v_out);
		s <<  "\nleg no. " << i+1 << ": " << m_seq[i-1]->get_name() << " to " << m_seq[i]->get_name() << std::endl; 
		s <<  "\tDuration: (days)" << T[i] << std::endl; 
		s <<  "\tFly-by epoch: " << t_P[i-1] << " (" << t_P[i-1].mjd2000() << " mjd2000) " << std::endl; 
		s <<  "\tFly-by altitude (km): " << (x[4*i+1]*m_seq[i-1]->get_radius()-m_seq[i-1]->get_radius())/1000.0 << std::endl; 
		s <<  "\tPlanet position (m): " << r_P[i-1] << std::endl; 
		s <<  "\tPlanet velocity (m/s): " << mem_vP << std::endl; 
		s <<  "\tV in (m/s): " << mem_vin << std::endl; 
		s <<  "\tV out (m/s): " << mem_vout << std::endl << std::endl;

		s <<  "\tDSM after (days): "  << x[4*i+2]*T[i] << std::endl; 
		s <<  "\tDSM magnitude (m/s): " << DV[i] << std::endl; 
		s <<  "\tClosest distance (JR): " << d << std::endl; 
		s <<  "\tApoapsis at closest distance (JR): " << ra << std::endl; 
 	}
	
	s << "\nArrival at " << m_seq[m_seq.size()-1]->get_name() << std::endl; 
	kep_toolbox::diff(v_out, v_end_l, v_P[m_seq.size()-1]);
	s <<  "Arrival epoch: "  << t_P[m_seq.size()-1] << " (" << t_P[m_seq.size()-1].mjd2000() << " mjd2000) " << std::endl; 
	s <<  "Arrival Vinf (m/s): " << v_out << " - " << kep_toolbox::norm(v_out) << std::endl; 
	s <<  "Total mission time (days): " << std::accumulate(T.begin(),T.end(),0.0) << std::endl; 
	return s.str();
}
std::string mga_incipit::get_name() const
{
	return "MGA-INCIPIT (CAPTURE AT JUPITER)";
}


/// Gets the times of flight
/**
 * @return[out] vector of times of flight 
 */
const std::vector<std::vector<double> >& mga_incipit::get_tof() const {
	return m_tof;
}

/// Sets the times of flight
/**
 * This setter changes the problem bounds as to define a minimum and a maximum allowed total time of flight
 *
 * @param[in] tof vector of times of flight 
 */
void mga_incipit::set_tof(const std::vector<std::vector<double> >& tof) {
	if (tof.size() != (m_seq.size())) {
		pagmo_throw(value_error,"The time-of-flight vector (tof) has the wrong length");  
	}
	m_tof = tof;
	for (size_t i=0; i< m_seq.size(); ++i) {
		set_bounds(3+4*i,tof[i][0],tof[i][1]);
	}
}

/// Gets the planetary sequence defining the interplanetary mga-1dsm mission
/**
 * @return An std::vector containing the kep_toolbox::planets
 */
std::vector<kep_toolbox::planet::planet_ptr> mga_incipit::get_sequence() const {
	return m_seq;
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the values vector, the weights vectors and the max weight. It is concatenated
 * with the base::problem human_readable
 */
std::string mga_incipit::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "\n\tSequence: ";
	for (size_t i = 0; i<m_seq.size() ;++i) {
		oss << m_seq[i]->get_name() << " ";
	}
	oss << "\n\tTime of flights?: ";
	for (size_t i=0; i<m_seq.size(); ++i) {
	  oss << m_tof[i]<<' ';
	}
	return oss.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::mga_incipit)
