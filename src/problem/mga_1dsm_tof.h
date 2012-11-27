/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#ifndef PAGMO_PROBLEM_MGA_1DSM_TOF_H
#define PAGMO_PROBLEM_MGA_1DSM_TOF_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "../types.h"
#include "base.h"
#include "../keplerian_toolbox/planet_ss.h"
#include "../keplerian_toolbox/epoch.h"


namespace pagmo{ namespace problem {

/// A generic MGA-1DSM Problem
/**
 *
 * This class defines the global optimization problem (box-bounded, continuous) of an interplanetary trajectory modelled
 * as a Multiple Gravity Assist mission allowing one only Deep Space Manouvre per leg.
 * 
 * The decision vector is [t0] + [u,v,Vinf,eta1,T1] + [beta, rp/rP, eta2,T2] ..... in the units: [mjd2000, days] + [nd,nd,km/s,nd,days] + [rad,nd,nd,days] + ....
 * where Vinf = Vinf_mag*(cos(theta)*cos(phi)i+cos(theta)*sin(phi)j+sin(phi)k) and theta = 2*pi*u and phi = acos(2*v-1)-pi/2
 *  
 * Each leg time-of-flight is directly encoded (as T1, T2, ...) in contrast to mga_1dsm_alpha. Thus you have to define the bounds on the time of 
 * flights for each leg separately upon construction.
 * 
 * NOTE: The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.
 *
 * @see Izzo: "Global Optimization and Space Pruning for Spacecraft Trajectory Design, Spacecraft Trajectory Optimization, Conway, B. (Eds.), Cambridge University Press, pp.178-199, 2010)
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE mga_1dsm_tof: public base
{
	public:
		mga_1dsm_tof(const std::vector<kep_toolbox::planet_ptr> = construct_default_sequence(), 
			 const kep_toolbox::epoch t0_l = kep_toolbox::epoch(0), const kep_toolbox::epoch t0_r = kep_toolbox::epoch(1000),
			 const std::vector<std::vector<double> > = construct_default_tof(),
			 const double vinf_l = 0.5, const double vinf_u = 2.5,
			 const bool mo = false, const bool add_vinf_dep = false, const bool add_vinf_arr = true);
		mga_1dsm_tof(const mga_1dsm_tof&);
		base_ptr clone() const;
		
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x) const;
		void set_tof(const std::vector<std::vector<double> >);
		void set_launch_window(const kep_toolbox::epoch&, const kep_toolbox::epoch&);
		void set_vinf(const double);
		std::vector<kep_toolbox::planet_ptr> get_sequence() const;
		std::vector<std::vector<double> > get_tof() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		std::string human_readable_extra() const;
		static const std::vector<kep_toolbox::planet_ptr> construct_default_sequence() {
			std::vector<kep_toolbox::planet_ptr> retval;
			retval.push_back(kep_toolbox::planet_ss("earth").clone());
			retval.push_back(kep_toolbox::planet_ss("venus").clone());
			retval.push_back(kep_toolbox::planet_ss("earth").clone());
			return retval;
		};
		static const std::vector<std::vector<double> > construct_default_tof() {
			std::vector<std::vector <double> > retval;
			std::vector<double> e2v;
			e2v.push_back(200.0);
			e2v.push_back(700.0);
			std::vector<double> v2e;
			e2v.push_back(200.0);
			e2v.push_back(700.0);
			retval.push_back(e2v);
			retval.push_back(v2e);
			return retval;
		};
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
		std::vector<kep_toolbox::planet_ptr> m_seq;
		const size_t m_n_legs;
		bool m_add_vinf_dep;
		bool m_add_vinf_arr;
};

}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::mga_1dsm_tof);
#endif // PAGMO_PROBLEM_MGA_1DSM_TOF_H
