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

#ifndef PAGMO_PROBLEM_MGA_INCIPIT_CSTRS_H
#define PAGMO_PROBLEM_MGA_INCIPIT_CSTRS_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"
#include <keplerian_toolbox/planet/gtoc6.h>
#include <keplerian_toolbox/epoch.h>


namespace pagmo{ namespace problem {

/// The beginning of the GTOC6 Jupiter Capture Trajectory
/**
 *
 * A PyGMO global optimization problem (box-bounded, continuous) representing a capture
 * in the Jupiter system. Two constraints are considered: 1. Closest distance, 2. Time of flight
 *
 * Decision vector:
 * [t0,u,v,T0] + [beta1, rp1/rP1, eta1,T1] + .... 
 * 
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE mga_incipit_cstrs: public base
{
	public:
		mga_incipit_cstrs(const std::vector<kep_toolbox::planet::planet_ptr> = construct_default_sequence(),
			 const kep_toolbox::epoch t0_l = kep_toolbox::epoch(7305.0),
			 const kep_toolbox::epoch t0_u = kep_toolbox::epoch(11323.0),
			 const std::vector<std::vector<double> > tof = construct_default_tofs(),
			 double Tmax = 365.25,
			 double Dmin = 0.2
			 );
		mga_incipit_cstrs(const mga_incipit_cstrs&);
		base_ptr clone() const;
		
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x) const;
		void set_tof(const std::vector<std::vector<double> >&);
		const std::vector<std::vector<double> >& get_tof() const;
		std::vector<kep_toolbox::planet::planet_ptr> get_sequence() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		std::string human_readable_extra() const;
		
	private:
		static const std::vector<kep_toolbox::planet::planet_ptr> construct_default_sequence() {
			std::vector<kep_toolbox::planet::planet_ptr> retval;
			retval.push_back(kep_toolbox::planet::gtoc6("io").clone());
			retval.push_back(kep_toolbox::planet::gtoc6("io").clone());
			retval.push_back(kep_toolbox::planet::gtoc6("europa").clone());
			return retval;
		}
		static const std::vector<std::vector<double> > construct_default_tofs() {
			std::vector<std::vector<double> > retval;
			std::vector<double> dumb(2);
			dumb[0] = 100;dumb[1] = 200;
			retval.push_back(dumb);
			dumb[0] = 3;dumb[1] = 200;
			retval.push_back(dumb);
			dumb[0] = 4;dumb[1] = 100;
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
			ar & m_tof;
			ar & const_cast<double &>(m_tmax);
			ar & const_cast<double &>(m_dmin);
		}
		std::vector<kep_toolbox::planet::planet_ptr> m_seq;
		std::vector<std::vector<double> > m_tof;
		const double m_tmax;
		const double m_dmin;
};

}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::mga_incipit_cstrs)
#endif // PAGMO_PROBLEM_MGA_INCIPIT_CSTRS_H
