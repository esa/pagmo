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

#ifndef PAGMO_PROBLEM_MGA_1DSM_TOF_H
#define PAGMO_PROBLEM_MGA_1DSM_TOF_H

#include <string>
#include <boost/array.hpp>
#include <keplerian_toolbox/planet/jpl_low_precision.h>
#include <keplerian_toolbox/epoch.h>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"



namespace pagmo{ namespace problem {

/// A generic MGA-1DSM Problem (tof encoding)
/**
 *
 * This class defines the global optimization problem (box-bounded, continuous) of an interplanetary trajectory modelled
 * as a Multiple Gravity Assist mission allowing one only Deep Space Manouvre per leg.
 * 
 * \image html mga_1dsm.gif "Visualization of an inetrplanetary trajectory to jupiter as encoded by the mga_1dsm"
 * \image latex mga_1dsm.png "Visualization of an inetrplanetary trajectory to jupiter as encoded by the mga_1dsm" width=5cm
 * 
 * The decision vector is \f$ [t_0] + [u,v,V_{\infty},\eta_1,T_1] + [\beta, r_p/r_P, \eta_2,T_2]\f$ ..... in the units: [mjd2000] + [nd,nd,km/s,nd,days] + [rad,nd,nd,days] + ....
 * where \f$ \mathbf V_{\infty} = V_{\infty}*(\cos(\theta)\cos(\phi)\mathbf i+\cos(\theta)\sin(\phi)\mathbf j+\sin(\phi)\mathbf k) \f$ and \f$ \theta = 2\pi u, \phi = acos(2v-1)-\pi/2 \f$
 *  
 * Each leg time-of-flight is directly encoded (as T1, T2, ...) in contrast to mga_1dsm_alpha. Thus you have to define the bounds on the time of 
 * flights for each leg separately upon construction.
 * 
 * NOTE: The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.
 *
 * @see Izzo: "Global Optimization and Space Pruning for Spacecraft Trajectory Design, Spacecraft Trajectory Optimization, Conway, B. (Eds.), Cambridge University Press, pp.178-199, 2010)
 * @author Dario Izzo (dario.izzo@esa.int)
 * @author Marcus Maertens (mmarcusx@gmail.com)
 */
class __PAGMO_VISIBLE mga_1dsm_tof: public base
{
	public:
		mga_1dsm_tof(const std::vector<kep_toolbox::planet::planet_ptr> = construct_default_sequence(), 
			 const kep_toolbox::epoch t0_l = kep_toolbox::epoch(0), const kep_toolbox::epoch t0_r = kep_toolbox::epoch(1000),
			 const std::vector<boost::array<double,2> > = construct_default_tof(),
			 const double vinf_l = 0.5, const double vinf_u = 2.5,
			 const bool mo = false, const bool add_vinf_dep = false, const bool add_vinf_arr = true);
		mga_1dsm_tof(const mga_1dsm_tof&);
		base_ptr clone() const;
		
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x, bool extended_output = false) const;
		void set_tof(const std::vector<boost::array<double,2> >);
		void set_launch_window(const kep_toolbox::epoch&, const kep_toolbox::epoch&);
		void set_vinf(const double);
		std::vector<kep_toolbox::planet::planet_ptr> get_sequence() const;
		std::vector<std::vector<double> > get_tof() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		std::string human_readable_extra() const;
		
	private:
		static const std::vector<kep_toolbox::planet::planet_ptr> construct_default_sequence() {
			std::vector<kep_toolbox::planet::planet_ptr> retval;
			retval.push_back(kep_toolbox::planet::jpl_lp("earth").clone());
			retval.push_back(kep_toolbox::planet::jpl_lp("venus").clone());
			retval.push_back(kep_toolbox::planet::jpl_lp("earth").clone());
			return retval;
		}
		static const std::vector<boost::array<double,2> > construct_default_tof() {
			std::vector<boost::array<double,2> > retval;
			boost::array<double,2> dumb = {{ 50,900 }};  
			retval.push_back(dumb);
			retval.push_back(dumb);
			return retval;
		}
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_seq;
			ar & const_cast<size_t &>(m_n_legs);
			ar & m_add_vinf_dep;
			ar & m_add_vinf_arr;
		}
		std::vector<kep_toolbox::planet::planet_ptr> m_seq;
		const size_t m_n_legs;
		bool m_add_vinf_dep;
		bool m_add_vinf_arr;
};

}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::mga_1dsm_tof)
#endif // PAGMO_PROBLEM_MGA_1DSM_TOF_H
