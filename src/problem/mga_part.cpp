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

#include "mga_part.h"


#define ASTRO_JR 71492000.0 //m

namespace pagmo { namespace problem {

 
/// Constructor
/**
 * Constructs a global optimization problem (box-bounded, continuous) representing a part og the gtoc6 jupiter moon tour
 *  
 * @param[in] seq std::vector of kep_toolbox::planet_ptr containing the encounter sequence for the trajectoty (including the initial planet)
 * @param[in] tof time-of-flight vector containing lower and upper bounds (in days) for the various legs time of flights
 * @param[in] t0 starting kep_toolbox::epoch
 * @param[in] vinf_in incoming velocity at the departing moon (before fly-by, absolute)
 * 
 * @throws value_error if the planets in seq do not all have the same central body gravitational constant
 * @throws value_error if tof has a size different from seq.size()
 * @throws value_error if seq.size() is < 2
 */
mga_part::mga_part(const std::vector<kep_toolbox::planet::planet_ptr> seq, 
			 const std::vector<std::vector<double> > tof,
			 const kep_toolbox::epoch t0,
			 const kep_toolbox::array3D vinf_in
			 ): 
			 base(4*(seq.size()-1)), m_tof(tof), m_t0(t0), m_vinf_in(vinf_in)
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
	if (tof.size() != (seq.size()-1)) {
		pagmo_throw(value_error,"The time-of-flight vector (tof) has the wrong length");  
	}
	for (size_t i = 0; i < tof.size(); ++i) {
		if (tof[i].size()!=2) pagmo_throw(value_error,"Each element of the time-of-flight vector (tof)  needs to have dimension 2 (lower and upper bound)"); 
	}
	//We check the sequence contains at least two planets
	if (seq.size() < 2) {
		pagmo_throw(value_error,"The sequence must contain at least 2 planets");  
	}
	
	// Filling in the planetary sequence data member. This requires to construct the polymorphic planets via their clone method 
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < seq.size(); ++i) {
		m_seq.push_back(seq[i]->clone());
	}
	
	// Now setting the problem bounds
	size_type dim(4*(seq.size()-1));
	decision_vector lb(dim), ub(dim);
	

	// all legs
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < (m_seq.size() - 1); ++i) {
		lb[4*i] = - 2 * boost::math::constants::pi<double>();    ub[4*i] = 2 * boost::math::constants::pi<double>();
		lb[4*i+1] = 1.1;  ub[4*i+1] = 30;
		lb[4*i+2] = 1e-5; ub[4*i+2] = 1-1e-5;
		lb[4*i+3] = m_tof[i][0]; ub[4*i+3] = m_tof[i][1];
	}
	
	// Adjusting the minimum and maximum allowed fly-by rp to the one defined in the kep_toolbox::planet class
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < (m_seq.size() - 1); ++i) {
		lb[4*i+1] = m_seq[i]->get_safe_radius() / m_seq[i]->get_radius();
		ub[4*i+1] = (m_seq[i]->get_radius() + 2000000) / m_seq[i]->get_radius(); //from gtoc6 problem description
	}
	set_bounds(lb,ub);
}

/// Copy Constructor. Performs a deep copy
mga_part::mga_part(const mga_part &p) : base(p.get_dimension()), m_tof(p.m_tof), m_t0(p.m_t0), m_vinf_in(p.m_vinf_in)
{
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < p.m_seq.size();++i) {
		m_seq.push_back(p.m_seq[i]->clone());
	}
	set_bounds(p.get_lb(),p.get_ub());
}

/// Clone method.
base_ptr mga_part::clone() const
{
	return base_ptr(new mga_part(*this));
}

/// Implementation of the objective function.
void mga_part::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
try {
	double d,ra,d2,ra2;
	double common_mu = m_seq[0]->get_mu_central_body();
	// 1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
	std::vector<double> T(m_seq.size()-1,0.0);
	
	for (size_t i = 0; i < (m_seq.size()-1); ++i) {
		T[i] = x[4*i+3];
	}
	// 2 - We compute the epochs and ephemerides of the planetary encounters
	std::vector<kep_toolbox::epoch>   t_P(m_seq.size());
	std::vector<kep_toolbox::array3D> r_P(m_seq.size());
	std::vector<kep_toolbox::array3D> v_P(m_seq.size());
	std::vector<double> DV(m_seq.size() - 1);
	for (size_t i = 0; i<m_seq.size(); ++i) {
		t_P[i] = kep_toolbox::epoch(m_t0.mjd2000() + std::accumulate(T.begin(),T.begin()+i,0.0));
		m_seq[i]->eph(t_P[i], r_P[i], v_P[i]);
	}

	// 3 - We loop over the legs
	kep_toolbox::array3D v_out,r,v;
	kep_toolbox::array3D v_beg_l,v_end_l;
	kep_toolbox::sum(v_end_l, m_vinf_in, v_P[0]);
	for (size_t i = 0; i<(m_seq.size()-1); ++i) {
		// Fly-by
		kep_toolbox::fb_prop(v_out, v_end_l, v_P[i], x[4*i+1] * m_seq[i]->get_radius(), x[4*i], m_seq[i]->get_mu_self());
		// s/c propagation before the DSM
		r = r_P[i];
		v = v_out;
		
		kep_toolbox::propagate_lagrangian(r,v,x[4*i+2]*T[i]*ASTRO_DAY2SEC,common_mu);
		kep_toolbox::closest_distance(d, ra, r_P[i], v_out, r, v, common_mu);
		
		// Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
		double dt = (1-x[4*i+2])*T[i]*ASTRO_DAY2SEC;
		kep_toolbox::lambert_problem l2(r,r_P[i+1],dt,common_mu,false,false);
		v_end_l = l2.get_v2()[0];
		v_beg_l = l2.get_v1()[0];
		kep_toolbox::closest_distance(d2,ra2,r,v_beg_l, r_P[i+1], v_end_l, common_mu);
		
		if (d < d2)
		{
			d = d/ASTRO_JR;
		} else {
			d = d2/ASTRO_JR;
		}

		// DSM occuring at time nu2*T2
		kep_toolbox::diff(v_out, v_beg_l, v);
		// Penalized if the JR is too small (<0.2)
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

std::string mga_part::pretty(const std::vector<double> &x) const {
  
	// We set the std output format
	std::ostringstream s;
	s.precision(15);
	s << std::scientific;
	
	double d,ra,d2,ra2;

	double common_mu = m_seq[0]->get_mu_central_body();
	// 1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
	std::vector<double> T(m_seq.size()-1,0.0);
	
	for (size_t i = 0; i<m_seq.size()-1; ++i) {
		T[i] = x[4*i+3];
	}
	// 2 - We compute the epochs and ephemerides of the planetary encounters
	std::vector<kep_toolbox::epoch>   t_P(m_seq.size());
	std::vector<kep_toolbox::array3D> r_P(m_seq.size());
	std::vector<kep_toolbox::array3D> v_P(m_seq.size());
	std::vector<double> DV(m_seq.size() - 1);
	for (size_t i = 0; i<m_seq.size(); ++i) {
		t_P[i] = kep_toolbox::epoch(m_t0.mjd2000() + std::accumulate(T.begin(),T.begin()+i,0.0));
		m_seq[i]->eph(t_P[i], r_P[i], v_P[i]);
	}

	// 4 - And we proceed with each successive leg (if any)
	kep_toolbox::array3D v_out,r,v,mem_vin,mem_vout,mem_vP;
	kep_toolbox::array3D v_beg_l,v_end_l;
	kep_toolbox::sum(v_end_l, m_vinf_in, v_P[0]);
	for (size_t i = 0; i<m_seq.size()-1; ++i) {
		// Fly-by
		kep_toolbox::fb_prop(v_out, v_end_l, v_P[i], x[4*i+1] * m_seq[i]->get_radius(), x[4*i], m_seq[i]->get_mu_self());
		// s/c propagation before the DSM
		r = r_P[i];
		v = v_out;
		
		mem_vout = v_out;
		mem_vin = v_end_l;
		mem_vP = v_P[i];
		
		kep_toolbox::propagate_lagrangian(r,v,x[4*i+2]*T[i]*ASTRO_DAY2SEC,common_mu);
		kep_toolbox::closest_distance(d, ra, r_P[i], v_out, r, v, common_mu);
		
		// Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
		double dt = (1-x[4*i+2])*T[i]*ASTRO_DAY2SEC;
		kep_toolbox::lambert_problem l2(r,r_P[i+1],dt,common_mu,false,false);
		v_end_l = l2.get_v2()[0];
		v_beg_l = l2.get_v1()[0];
		kep_toolbox::closest_distance(d2,ra2,r,v_beg_l, r_P[i+1], v_end_l, common_mu);
		
		if (d < d2)
		{
			d = d/ASTRO_JR;
			ra = ra/ASTRO_JR;
		} else {
			d = d2/ASTRO_JR;
			ra = ra2/ASTRO_JR;
		}
		// DSM occuring at time nu2*T2
		// DSM occuring at time nu2*T2
		kep_toolbox::diff(v_out, v_beg_l, v);
		DV[i] = kep_toolbox::norm(v_out);
		s <<  "\nleg no. " << i+1 << ": " << m_seq[i]->get_name() << " to " << m_seq[i+1]->get_name() << std::endl; 
		s <<  "\tDuration: (days)" << T[i] << std::endl; 
		s <<  "\tFly-by epoch: " << t_P[i] << " (" << t_P[i].mjd2000() << " mjd2000) " << std::endl; 
		s <<  "\tFly-by altitude (km): " << (x[4*i+1]*m_seq[i]->get_radius()-m_seq[i]->get_radius())/1000.0 << std::endl; 
		s <<  "\tPlanet position (m): " << r_P[i] << std::endl; 
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
std::string mga_part::get_name() const
{
	return "MGA-PART (a part of the jupiter moon tour)";
}

/// Sets the times of flight
/**
 * This setter changes the problem bounds as to define a minimum and a maximum allowed total time of flight
 *
 * @param[in] tof vector of times of flight 
 */
void mga_part::set_tof(const std::vector<std::vector<double> >& tof) {
	if (tof.size() != (m_seq.size()-1)) {
		pagmo_throw(value_error,"The time-of-flight vector (tof) has the wrong length");  
	}
	m_tof = tof;
	for (size_t i=0; i< (m_seq.size()-1); ++i) {
		set_bounds(3+4*i,tof[i][0],tof[i][1]);
	}
}

/// Gets the times of flight
/**
 * @return const reference to m_tof
 */
const std::vector<std::vector<double> >& mga_part::get_tof() const {
	return m_tof;
}

/// Sets the betas
/**
 * This setter changes the problem bounds as to define a minimum and a maximum allowed beta for
 * each leg. Remember that beta controls the inclination of the planetocentric hyperbola
 *
 * @param[in] betas vector of bounds for betas (rad)
 */
void mga_part::set_betas(const std::vector<std::vector<double> >& betas) {
	if (betas.size() != (m_seq.size()-1)) {
		pagmo_throw(value_error,"The betas vector (betas) has the wrong length");  
	}
	for (size_t i=0; i< (m_seq.size()-1); ++i) {
		if (betas[i].size() !=2) {
			pagmo_throw(value_error,"Lower and upper bound vector needs to have dimension 2");  
		}
		set_bounds(4*i,betas[i][0],betas[i][1]);
	}
}

/// Gets the betas
/**
 * @return vector of betas bounds
 */
std::vector<std::vector<double> > mga_part::get_betas() const {
    std::vector<std::vector<double> > retval;
	std::vector<double> tmp(2);
	for (size_t i=0; i< (m_seq.size()-1); ++i) {
		tmp[0] = get_lb()[4*i];
		tmp[1] = get_ub()[4*i];
		retval.push_back(tmp);
	}
	return retval;
}

/// Sets the peri-planet bounds
/**
 * This setter changes the problem bounds as to define a minimum and a maximum allowed periplanet distance
 *
 * @param[in] rps vector of peri-planets (altitudes in km)
 */
void mga_part::set_rps(const std::vector<std::vector<double> >& rps) {
	if (rps.size() != (m_seq.size()-1)) {
		pagmo_throw(value_error,"The periplanets vector (rps) has the wrong length");  
	}
	for (size_t i=0; i< (m_seq.size()-1); ++i) {
		if (rps[i].size() !=2) {
			pagmo_throw(value_error,"Lower and upper bound vector needs to have dimension 2");  
		}
		set_bounds(1+4*i,(rps[i][0]*1000+m_seq[i]->get_radius())/m_seq[i]->get_radius(),(rps[i][1]*1000+m_seq[i]->get_radius())/m_seq[i]->get_radius());
	}
}

/// Gets the peri-planet bounds
/**
 * @return vector of periplanets (altitudes in km)
 */
std::vector<std::vector<double> > mga_part::get_rps() const {
	std::vector<std::vector<double> > retval;
	std::vector<double> tmp(2);
	for (size_t i=0; i< (m_seq.size()-1); ++i) {
		tmp[0] = ((get_lb()[4*i+1]*m_seq[i]->get_radius())-m_seq[i]->get_radius())/1000.0;
		tmp[1] = ((get_ub()[4*i+1]*m_seq[i]->get_radius())-m_seq[i]->get_radius())/1000.0;
		retval.push_back(tmp);
	}
	return retval;
}

/// Sets the start epoch
/**
 * This setter changes the starting epoch
 *
 * @param[in] t0 new start epoch
 */
void mga_part::set_t0(const kep_toolbox::epoch& t0) {
	m_t0 = t0;
}

/// Gets the start epoch
/**
 * @return start epoch
 */
const kep_toolbox::epoch& mga_part::get_t0() const {
	return m_t0;
}

/// Sets the start velocity 
/**
 * This setter changes the relative velocity of the spacecraft used to define the starting conditions on the traj
 *
 * @param[in] vinf_in new incoming relative velocity (m/s)
 */
void mga_part::set_vinf_in(const kep_toolbox::array3D& vinf_in) {
	m_vinf_in = vinf_in;
}

/// Gets the start velocity 
/**
 * @return relative velocity (m/s)
 */
const kep_toolbox::array3D& mga_part::get_vinf_in() const {
	return m_vinf_in;
}


/// Gets the planetary sequence defining the interplanetary mga-1dsm mission
/**
 * @return An std::vector containing the kep_toolbox::planets
 */
std::vector<kep_toolbox::planet::planet_ptr> mga_part::get_sequence() const {
	return m_seq;
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the values vector, the weights vectors and the max weight. It is concatenated
 * with the base::problem human_readable
 */
std::string mga_part::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "\n\tSequence: ";
	for (size_t i = 0; i<m_seq.size() ;++i) {
		oss << m_seq[i]->get_name() << " ";
	}
	oss << "\n\tTime of flights: ";
	for (size_t i=0; i<(m_seq.size()-1); ++i) {
	  oss << m_tof[i]<<' ';
	}
	oss << "\n\tStarting epoch: " << m_t0;
	oss << "\n\tStarting Vinf in: " << m_vinf_in << std::endl;
	return oss.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::mga_part)
