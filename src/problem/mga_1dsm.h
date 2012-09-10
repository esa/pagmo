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

#ifndef PAGMO_PROBLEM_MGA_1DSM_H
#define PAGMO_PROBLEM_MGA_1DSM_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "../types.h"
#include "base.h"
#include "../keplerian_toolbox/keplerian_toolbox.h"


namespace pagmo{ namespace problem {

/// A generic MGA-1DSM Problem
/**
 *
 * This class represents a global optimization problem (box-bounded, continuous) relative to an interplanetary trajectory modelled
 * as a Multiple Gravity Assist trajectory that allows one only Deep Space Manouvre between each leg.
 * 
 *
 * The decision vector is [t0,u,v,Vinf,eta1,T] + [beta, rp/rV, eta2,a2] ..... in the units: [mjd2000,nd,nd,km/s,nd,years] + [rad,nd,nd,nd] + ....
 * where Vinf = Vinf_mag*(cos(theta)*cos(phi)i+cos(theta)*sin(phi)j+sin(phi)k) and theta = 2*pi*u and phi = acos(2*v-1)-pi/2
 * Each leg time-of-flight can be obtained as Tn = T*an, T(n-1) = (T - Tn)*a(n-1), .... , Ti = (T-T(i+1)-T(i+2)- .... - Tn)*ai
 * 
 * NOTE: The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.
 *
 *
 * @see Izzo: "Global Optimization and Space Pruning for Spacecraft Trajectory Design, Spacecraft Trajectory Optimization, Conway, B. (Eds.), Cambridge University Press, pp.178-199, 2010)
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE mga_1dsm: public base
{
	public:
		mga_1dsm(const std::vector<kep_toolbox::planet_ptr> = construct_default_sequence(), 
			 const kep_toolbox::epoch t0_l = kep_toolbox::epoch(0), const kep_toolbox::epoch t0_r = kep_toolbox::epoch(1000),
			 const double tof_l = 1.0, const double tof_u = 5.0, 
			 const double vinf_l = 0.5, const double vinf_u = 2.5, 
			 const bool mo = false, const bool add_vinf = false);
		mga_1dsm(const mga_1dsm&);
		base_ptr clone() const;
		
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x) const;
		void set_tof(const double, const double);
		void set_launch_window(const kep_toolbox::epoch&, const kep_toolbox::epoch&);
		void set_vinf(const double);
		std::vector<kep_toolbox::planet_ptr> get_sequence() const;
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
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_seq;
			ar & const_cast<size_t &>(m_n_legs);
		}
		std::vector<kep_toolbox::planet_ptr> m_seq;
		const size_t m_n_legs;
		bool m_add_vinf;
};

}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::mga_1dsm);
#endif // PAGMO_PROBLEM_MGA_1DSM_H
