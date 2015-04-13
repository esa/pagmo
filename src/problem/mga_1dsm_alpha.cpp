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
#include <keplerian_toolbox/lambert_problem.h>

#include "mga_1dsm_alpha.h"


namespace pagmo { namespace problem {

 
/// Constructor
/**
 * Constructs a global optimization problem (box-bounded, continuous) representing an interplanetary trajectory modelled
 * as a Multiple Gravity Assist trajectory that allows one only Deep Space Manouvre between each leg.
 *  
 * @param[in] seq std::vector of kep_toolbox::planet_ptr containing the encounter sequence for the trajectoty (including the initial planet)
 * @param[in] t0_l kep_toolbox::epoch representing the lower bound for the launch epoch
 * @param[in] t0_u kep_toolbox::epoch representing the upper bound for the launch epoch
 * @param[in] tof_l lower bound for the mission duration (in days)
 * @param[in] tof_u upper bound for the mission duration (in days)
 * @param[in] vinf_l  minimum launch hyperbolic velocity allowed (in km/s)
 * @param[in] vinf_u  maximum launch hyperbolic velocity allowed (in km/s)
 * @param[in] mo: when true defines the problem as a multi-objective problem, returning total DV and time of flight
 * @param[in] add_vinf_dep: when true the objective functions contains also the contribution of the launch hypebolic velocity
 * @param[in] add_vinf_arr: when true the objective functions contains also the contribution of the arrival hypebolic velocity
 * 
 * @throws value_error if the planets in seq do not all have the same central body gravitational constant
 */
mga_1dsm_alpha::mga_1dsm_alpha(const std::vector<kep_toolbox::planet::planet_ptr> seq, 
			 const kep_toolbox::epoch t0_l, const kep_toolbox::epoch t0_u,
			 const double tof_l, const double tof_u, 
			 const double vinf_l, const double vinf_u, 
			 const bool mo, const bool add_vinf_dep, const bool add_vinf_arr) : 
			 base( 7 +  (int)(seq.size()-2) * 4, 0, 1 + (int)mo,0,0,0.0), m_seq(), m_n_legs(seq.size()-1), m_add_vinf_dep(add_vinf_dep), m_add_vinf_arr(add_vinf_arr)
{
	// We check that all planets have equal central body
	std::vector<double> mus(seq.size());
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i< seq.size(); ++i) {
		mus[i] = seq[i]->get_mu_central_body();
	}
	if ((unsigned int)std::count(mus.begin(), mus.end(), mus[0]) != mus.size()) {
		pagmo_throw(value_error,"The planets do not all have the same mu_central_body");  
	}
	// Filling in the planetary sequence data member. This requires to construct the polymorphic planets via their clone method 
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < seq.size(); ++i) {
		m_seq.push_back(seq[i]->clone());
	}
	
	// Now setting the problem bounds
	size_type dim(7 +  (m_n_legs-1) * 4);
	decision_vector lb(dim), ub(dim);
	
	// First leg
	lb[0] = t0_l.mjd2000(); ub[0] = t0_u.mjd2000();
	lb[1] = tof_l; ub[1] = tof_u;
	lb[2] = 0; lb[3] = 0; ub[2] = 1; ub[3] = 1;
	lb[4] = vinf_l * 1000; ub[4] = vinf_u * 1000;
	lb[5] = 1e-5; ub[5] = 1-1e-5;
	lb[6] = 1e-5; ub[6] = 1-1e-5;
	
	// Successive legs
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < m_n_legs - 1; ++i) {
		lb[7+4*i] = - 2 * boost::math::constants::pi<double>();    ub[7+4*i] = 2 * boost::math::constants::pi<double>();
		lb[8+4*i] = 1.1;  ub[8+4*i] = 100;
		lb[9+4*i] = 1e-5; ub[9+4*i] = 1-1e-5;
		lb[10+4*i] = 1e-5; ub[10+4*i] = 1-1e-5;
	}
	
	// Adjusting the minimum allowed fly-by rp to the one defined in the kep_toolbox::planet class
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 1; i < m_n_legs; ++i) {
		lb[4 + 4*i] = m_seq[i]->get_safe_radius() / m_seq[i]->get_radius();
	}
	set_bounds(lb,ub);
}

/// Copy Constructor. Performs a deep copy
mga_1dsm_alpha::mga_1dsm_alpha(const mga_1dsm_alpha &p) : base(p.get_dimension(), 0, p.get_f_dimension(),0,0,0.0), m_seq(), m_n_legs(p.m_n_legs), m_add_vinf_dep(p.m_add_vinf_dep), m_add_vinf_arr(p.m_add_vinf_arr) 
{
	for (std::vector<kep_toolbox::planet::planet_ptr>::size_type i = 0; i < p.m_seq.size();++i) {
		m_seq.push_back(p.m_seq[i]->clone());
	}
	set_bounds(p.get_lb(),p.get_ub());
}

/// Clone method.
base_ptr mga_1dsm_alpha::clone() const
{
	return base_ptr(new mga_1dsm_alpha(*this));
}

/// Implementation of the objective function.
void mga_1dsm_alpha::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
try {
	double common_mu = m_seq[0]->get_mu_central_body();
	// 1 - we 'decode' the chromosome recording the various times of flight (days) in the list T
	std::vector<double> T(m_n_legs,0.0);
	double alpha_sum = 0;
	
	for (size_t i = 0; i < m_n_legs; ++i) {
		double tmp = -log(x[6+4*i]);
		alpha_sum+= tmp;
		T[i] = x[1] * tmp;
	}
	for (size_t i = 0; i < m_n_legs; ++i) {
		T[i] /= alpha_sum;
	}

	// 2 - We compute the epochs and ephemerides of the planetary encounters
	std::vector<kep_toolbox::epoch>   t_P(m_n_legs + 1);
	std::vector<kep_toolbox::array3D> r_P(m_n_legs + 1);
	std::vector<kep_toolbox::array3D> v_P(m_n_legs + 1);
	std::vector<double> DV(m_n_legs + 1);
	for (size_t i = 0; i<(m_n_legs + 1); ++i) {
		t_P[i] = kep_toolbox::epoch(x[0] + std::accumulate(T.begin(), T.begin()+i, 0.0));
		m_seq[i]->eph(t_P[i], r_P[i], v_P[i]);
	}

	// 3 - We start with the first leg
	double theta = 2*boost::math::constants::pi<double>()*x[2];
	double phi = acos(2*x[3]-1)-boost::math::constants::pi<double>() / 2;

	kep_toolbox::array3D Vinf = { {x[4]*cos(phi)*cos(theta), x[4]*cos(phi)*sin(theta), x[4]*sin(phi)} };
	kep_toolbox::array3D v0;
	kep_toolbox::sum(v0, v_P[0], Vinf);
	kep_toolbox::array3D r(r_P[0]), v(v0);
	kep_toolbox::propagate_lagrangian(r,v,x[5]*T[0]*ASTRO_DAY2SEC,common_mu);

	// Lambert arc to reach seq[1]
	double dt = (1-x[5])*T[0]*ASTRO_DAY2SEC;
	kep_toolbox::lambert_problem l(r,r_P[1],dt,common_mu);
	kep_toolbox::array3D v_end_l = l.get_v2()[0];
	kep_toolbox::array3D v_beg_l = l.get_v1()[0];

	// First DSM occuring at time nu1*T1
	kep_toolbox::diff(v, v_beg_l, v);
	DV[0] = kep_toolbox::norm(v);

	// 4 - And we proceed with each successive leg (if any)
	kep_toolbox::array3D v_out;
	for (size_t i = 1; i<m_n_legs; ++i) {
		// Fly-by
		kep_toolbox::fb_prop(v_out, v_end_l, v_P[i], x[8+(i-1)*4] * m_seq[i]->get_radius(), x[7+(i-1)*4], m_seq[i]->get_mu_self());
		// s/c propagation before the DSM
		r = r_P[i];
		v = v_out;
		kep_toolbox::propagate_lagrangian(r,v,x[9+(i-1)*4]*T[i]*ASTRO_DAY2SEC,common_mu);

		// Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
		dt = (1-x[9+(i-1)*4])*T[i]*ASTRO_DAY2SEC;
		kep_toolbox::lambert_problem l2(r,r_P[i+1],dt,common_mu);
	  	v_end_l = l2.get_v2()[0];
		v_beg_l = l2.get_v1()[0];

		// DSM occuring at time nu2*T2
		kep_toolbox::diff(v, v_beg_l, v);
		DV[i] = kep_toolbox::norm(v);
	}
	
	// Last Delta-v
	kep_toolbox::diff(v, v_end_l, v_P[m_n_legs]);
	DV[m_n_legs] = kep_toolbox::norm(v);
	
	
	// Now we return the objective(s) function
	f[0] = std::accumulate(DV.begin(),DV.end()-1,0.0); 
	if (m_add_vinf_dep) {
		f[0] += x[4];
	}
	if (m_add_vinf_arr) {
		f[0] += DV[DV.size()-1];
	}
	if (get_f_dimension() == 2){
		f[1] = std::accumulate(T.begin(),T.end(),0.0);
	} 
//Here the lambert solver or the lagrangian propagator went wrong
} catch (...) {
	f[0] = boost::numeric::bounds<double>::highest();
	if (get_f_dimension() == 2){
		f[1] = boost::numeric::bounds<double>::highest();
	}
} 
}

/// Outputs a stream with the trajectory data
/**
 * While the chromosome contains all necessary information to describe a trajectory, mission analysis
 * often require a different set of data to evaluate its use. This method outputs a stream with
 * information on the trajectory that is otherwise 'hidden' in the chromosome
 *
 * \param[in] x chromosome representing the trajectory in the optimization process
 * \returns an std::string with launch dates, DV magnitues and other information on the trajectory
 */


std::string mga_1dsm_alpha::pretty(const std::vector<double> &x) const {
  
	double common_mu = m_seq[0]->get_mu_central_body();
	// We set the std output format
	std::ostringstream s;
	s.precision(15);
	s << std::scientific;

	// 1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
	std::vector<double> T(m_n_legs,0.0);
	double alpha_sum = 0;
	
	for (size_t i = 0; i < m_n_legs; ++i) {
		double tmp = -log(x[6+4*i]);
		alpha_sum+= tmp;
		T[i] = x[1] * tmp;
	}
	for (size_t i = 0; i < m_n_legs; ++i) {
		T[i] /= alpha_sum;
	}

	// 2 - We compute the epochs and ephemerides of the planetary encounters
	std::vector<kep_toolbox::epoch>   t_P(m_n_legs + 1);
	std::vector<kep_toolbox::array3D> r_P(m_n_legs + 1);
	std::vector<kep_toolbox::array3D> v_P(m_n_legs + 1);
	std::vector<double> DV(m_n_legs + 1);
	for (size_t i = 0; i<(m_n_legs + 1); ++i) {
		t_P[i] = kep_toolbox::epoch(x[0] + std::accumulate(T.begin(), T.begin()+i, 0.0));
		m_seq[i]->eph(t_P[i], r_P[i], v_P[i]);
	}

	// 3 - We start with the first leg
	s << "First Leg: " <<  m_seq[0]->get_name() << " to " << m_seq[1]->get_name() << std::endl; 
	s << "Departure: " << t_P[0] << " (" << t_P[0].mjd2000() << " mjd2000)" << std::endl;
	s << "Duration: " << T[0] << " days" << std::endl;
	s << "VINF: " << x[4] / 1000 << " km/sec" << std::endl;
	s << "DSM after " << x[5]*T[0] << " days" << std::endl;
	double theta = 2*boost::math::constants::pi<double>()*x[2];
	double phi = acos(2*x[3]-1)-boost::math::constants::pi<double>() / 2;

	kep_toolbox::array3D Vinf = { {x[4]*cos(phi)*cos(theta), x[4]*cos(phi)*sin(theta), x[4]*sin(phi)} };

	kep_toolbox::array3D v0;
	kep_toolbox::sum(v0, v_P[0], Vinf);
	kep_toolbox::array3D r(r_P[0]), v(v0);
	kep_toolbox::propagate_lagrangian(r,v,x[5]*T[0]*ASTRO_DAY2SEC,common_mu);

	// Lambert arc to reach seq[1]
	double dt = (1-x[5])*T[0]*ASTRO_DAY2SEC;
	kep_toolbox::lambert_problem l(r,r_P[1],dt,common_mu);
	kep_toolbox::array3D v_end_l = l.get_v2()[0];
	kep_toolbox::array3D v_beg_l = l.get_v1()[0];

	// First DSM occuring at time nu1*T1
	kep_toolbox::diff(v, v_beg_l, v);
	DV[0] = kep_toolbox::norm(v);
	s << "DSM magnitude: " << DV[0] << " m/s"  << std::endl;

	// 4 - And we proceed with each successive leg (if any)
	kep_toolbox::array3D v_out;
	for (size_t i = 1; i<m_n_legs; ++i) {
		s << "\nleg no: " << i+1 << ": " << m_seq[i]->get_name() << " to " << m_seq[i+1]->get_name()  << std::endl;
		s << "Duration: " << T[i] << " days"  << std::endl;
		s << "Fly-by epoch: " << t_P[i] << " (" << t_P[i].mjd2000() << " mjd2000) "  << std::endl;
		s << "Fly-by radius: " << x[8+(i-1)*4] << " planetary radii" << std::endl;
		s << "DSM after " << x[9+(i-1)*4]*T[i] << " days" << std::endl;
		// Fly-by
		kep_toolbox::fb_prop(v_out, v_end_l, v_P[i], x[8+(i-1)*4] * m_seq[i]->get_radius(), x[7+(i-1)*4], m_seq[i]->get_mu_self());
		// s/c propagation before the DSM
		r = r_P[i];
		v = v_out;
		kep_toolbox::propagate_lagrangian(r,v,x[9+(i-1)*4]*T[i]*ASTRO_DAY2SEC,common_mu);

		// Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
		dt = (1-x[9+(i-1)*4])*T[i]*ASTRO_DAY2SEC;
		kep_toolbox::lambert_problem l2(r,r_P[i+1],dt,common_mu);
	  	v_end_l = l2.get_v2()[0];
		v_beg_l = l2.get_v1()[0];

		// DSM occuring at time nu2*T2
		kep_toolbox::diff(v, v_beg_l, v);
		DV[i] = kep_toolbox::norm(v);
		s << "DSM magnitude: " << DV[i] << "m/s" << std::endl;
	}
	
	// Last Delta-v
	kep_toolbox::diff(v, v_end_l, v_P[m_n_legs]);
	DV[m_n_legs] = kep_toolbox::norm(v);
	s << "\nArrival at " << m_seq[m_n_legs]->get_name() << std::endl;
	s << "Arrival epoch: " << t_P[m_n_legs] << " (" << t_P[m_n_legs].mjd2000() << " mjd2000) " << std::endl;
	s << "Arrival Vinf: " << DV[m_n_legs] << "m/s" << std::endl;
	s << "Total mission time: " << std::accumulate(T.begin(),T.end(),0.0)/365.25 << " years" << std::endl;
	return s.str();
}
std::string mga_1dsm_alpha::get_name() const
{
	return "MGA-1DSM (alpha-encoding)";
}

/// Sets the mission time of flight
/**
 * This setter changes the problem bounds as to define a minimum and a maximum allowed total time of flight
 *
 * @param[in] tl Lower bownd on the time of flight in days
 * @param[in] tu Upper bownd on the time of flight in days
 */
void mga_1dsm_alpha::set_tof(double tl, double tu) {
	set_bounds(1,tl,tu);
}

/// Sets the mission launch window
/**
 * This setter changes the problem bounds as to define the launch window
 *
 * @param[in] start starting epoch of the launch window
 * @param[in] end final epoch of the launch window
 */
void mga_1dsm_alpha::set_launch_window(const kep_toolbox::epoch& start, const kep_toolbox::epoch& end) {
	set_bounds(0,start.mjd2000(),end.mjd2000());
}

/// Sets the launch hyperbolic velocity
/**
 * This setter changes the problem bounds as to define the maximum allowed launch hyperbolic velocity
 *
 * @param[in] vinf maximum allowed hyperbolic velocity in km/sec
 */
void mga_1dsm_alpha::set_vinf(const double vinf) {
	set_ub(4,vinf*1000);
}

/// Gets the sequence of time of flights for the mga-1dsm mission
/**
 * @return An std::vector of two doubles containing the lower and upper bound on total time of flight
 */
std::vector<double> mga_1dsm_alpha::get_tof() const {
	std::vector<double> retval;
	const decision_vector lb = get_lb();
	const decision_vector ub = get_ub();
	retval.push_back(lb[1]);
	retval.push_back(ub[1]);
	return retval;
}


/// Gets the planetary sequence defining the interplanetary mga-1dsm mission
/**
 * @return An std::vector containing the kep_toolbox::planets
 */
std::vector<kep_toolbox::planet::planet_ptr> mga_1dsm_alpha::get_sequence() const {
	return m_seq;
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the values vector, the weights vectors and the max weight. It is concatenated
 * with the base::problem human_readable
 */
std::string mga_1dsm_alpha::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "\n\tSequence: ";
	for (size_t i = 0; i<m_seq.size() ;++i) {
		oss << m_seq[i]->get_name() << " ";
	}
	oss << "\n\tAdd launcher vinf to the objective?: " << (m_add_vinf_dep? " True":" False") << std::endl;
	oss << "\n\tAdd arrival vinf to the objective?: " << (m_add_vinf_arr? " True":" False") << std::endl;
	return oss.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::mga_1dsm_alpha)
