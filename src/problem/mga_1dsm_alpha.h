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

#ifndef PAGMO_PROBLEM_MGA_1DSM_ALPHA_H
#define PAGMO_PROBLEM_MGA_1DSM_ALPHA_H

#include <string>
#include <keplerian_toolbox/planet/jpl_low_precision.h>
#include <keplerian_toolbox/epoch.h>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"



namespace pagmo{ namespace problem {

/// A generic MGA-1DSM Problem
/**
 *
 * This class defines the global optimization problem (box-bounded, continuous) of an interplanetary trajectory modelled
 * as a Multiple Gravity Assist mission allowing one only Deep Space Manouvre per leg.
 * 
 * \image html mga_1dsm.gif "Visualization of an inetrplanetary trajectory to jupiter as encoded by the mga_1dsm"
 * \image latex mga_1dsm.png "Visualization of an inetrplanetary trajectory to jupiter as encoded by the mga_1dsm" width=5cm
 * 
 * The decision vector is \f$ [t_0,T] + [u,v,V_{\infty},\eta_1,\alpha_1] + [\beta, r_p/r_P, \eta_2,\alpha_2]\f$ ..... in the units: [mjd2000, days] + [nd,nd,km/s,nd,nd] + [rad,nd,nd,nd] + ....
 * where \f$ \mathbf V_{\infty} = V_{\infty}*(\cos(\theta)\cos(\phi)\mathbf i+\cos(\theta)\sin(\phi)\mathbf j+\sin(\phi)\mathbf k) \f$ and \f$ \theta = 2\pi u, \phi = acos(2v-1)-\pi/2 \f$
 * 
 * Each leg time-of-flight can be obtained as \f$ T_n = T log(\alpha_n) / \sum_i(log(\alpha_i))\f$. This is
 * what we call \f$\alpha\f$-encoding as opposed to the tof encoding implemented in mga_1dsm_tof
 * 
 * This encoding allows the optimizer more flexibility in choosing the flybys and creates a
 * better multiobjective problem (the total tof can be set), but might greater a more difficult problem. 
 * The probability of a leg having a duration smaller or larger than some T, is the same for each planet in the sequence.
 * 
 * NOTE: The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.
 *
 * @see Izzo: "Global Optimization and Space Pruning for Spacecraft Trajectory Design, Spacecraft Trajectory Optimization, Conway, B. (Eds.), Cambridge University Press, pp.178-199, 2010)
 * @author Dario Izzo (dario.izzo@esa.int)
 * @author Marcus Maertens (mmarcusx@gmail.com)
 */
class __PAGMO_VISIBLE mga_1dsm_alpha: public base
{
	public:
		mga_1dsm_alpha(const std::vector<kep_toolbox::planet::planet_ptr> = construct_default_sequence(), 
			 const kep_toolbox::epoch t0_l = kep_toolbox::epoch(0), const kep_toolbox::epoch t0_r = kep_toolbox::epoch(1000),
			 const double tof_l = 1.0*365.25, const double tof_u = 5.0*365.25, 
			 const double vinf_l = 0.5, const double vinf_u = 2.5, 
			 const bool mo = false, const bool add_vinf_dep = false, const bool add_vinf_arr = true);
		mga_1dsm_alpha(const mga_1dsm_alpha&);
		base_ptr clone() const;
		
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x) const;
		void set_tof(const double, const double);
		void set_launch_window(const kep_toolbox::epoch&, const kep_toolbox::epoch&);
		void set_vinf(const double);
		std::vector<kep_toolbox::planet::planet_ptr> get_sequence() const;
		std::vector<double> get_tof() const;
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

BOOST_CLASS_EXPORT_KEY(pagmo::problem::mga_1dsm_alpha)
#endif // PAGMO_PROBLEM_MGA_1DSM_ALPHA_H
