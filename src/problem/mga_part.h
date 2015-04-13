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

#ifndef PAGMO_PROBLEM_MGA_PART_H
#define PAGMO_PROBLEM_MGA_PART_H

#include <string>
#include <keplerian_toolbox/planet/gtoc6.h>
#include <keplerian_toolbox/astro_constants.h>
#include <keplerian_toolbox/epoch.h>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"



namespace pagmo{ namespace problem {

/// A part of the GTOC6 Jupiter Capture Trajectory
/**
 *
 * A PyGMO global optimization problem (box-bounded, continuous) representing a part of the gtoc6 preliminary trajectory design
	
 * Decision vector:
 * [beta1, rp1/rP1, eta1,T1] + [beta2, rp2/rP2, eta2,T2] + ....
 * 
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE mga_part: public base
{
	public:
		mga_part(const std::vector<kep_toolbox::planet::planet_ptr> = construct_default_sequence(), 
			 const std::vector<std::vector<double> > tof = construct_default_tofs(),
			 const kep_toolbox::epoch t0 = kep_toolbox::epoch(11000),
			 const kep_toolbox::array3D vinf_in = construct_default_v()
			 );
		mga_part(const mga_part&);
		base_ptr clone() const;
		
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x) const;
		
		void set_tof(const std::vector<std::vector<double> >&);
		const std::vector<std::vector<double> >& get_tof() const;
		void set_t0(const kep_toolbox::epoch&);
		const kep_toolbox::epoch& get_t0() const;
		void set_vinf_in(const kep_toolbox::array3D&);
		const kep_toolbox::array3D& get_vinf_in() const;
		void set_betas(const std::vector<std::vector<double> >&);
		std::vector<std::vector<double> > get_betas() const;
		void set_rps(const std::vector<std::vector<double> >&);
		std::vector<std::vector<double> > get_rps() const;
		std::vector<kep_toolbox::planet::planet_ptr> get_sequence() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		std::string human_readable_extra() const;
	private:
		static const std::vector<kep_toolbox::planet::planet_ptr> construct_default_sequence() {
			std::vector<kep_toolbox::planet::planet_ptr> retval;
			retval.push_back(kep_toolbox::planet::gtoc6("europa").clone());
			retval.push_back(kep_toolbox::planet::gtoc6("europa").clone());
			retval.push_back(kep_toolbox::planet::gtoc6("europa").clone());
			return retval;
		}
		static const kep_toolbox::array3D construct_default_v() {
			const kep_toolbox::array3D retval = { {1500.0,2350.0,145.0} };
			return retval;
		}
		static const std::vector<std::vector<double> > construct_default_tofs() {
			std::vector<std::vector<double> > retval;
			std::vector<double> dumb(2);
			dumb[0] = 10;dumb[1] = 40;
			retval.push_back(dumb);
			dumb[0] = 10;dumb[1] = 40;
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
			ar & m_t0;
			ar & m_vinf_in;
		}
		std::vector<kep_toolbox::planet::planet_ptr> 		m_seq;
		std::vector<std::vector<double> > 					m_tof;
		kep_toolbox::epoch									m_t0;
		kep_toolbox::array3D								m_vinf_in;
};

}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::mga_part)
#endif // PAGMO_PROBLEM_MGA_PART_H
