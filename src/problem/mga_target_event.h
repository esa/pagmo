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

#ifndef PAGMO_PROBLEM_MGA_TARGET_EVENT_H
#define PAGMO_PROBLEM_MGA_TARGET_EVENT_H

#include <string>
#include <keplerian_toolbox/epoch.h>
#include <keplerian_toolbox/planet/base.h>
#include <keplerian_toolbox/planet/jpl_low_precision.h>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"



namespace pagmo{ namespace problem {

/// Used in GTOC 7
/**
 *
 * A PyGMO global optimization problem (box-bounded, continuous)
 * representing the problem of targeting a given event in the
 * overall mission timeline. An event is defined by a pagmo::planet
 * and an epoch.

 * Decision vector:
 *      [T,tw,Vinf,u,v,eta] in [days, n/a, m/s, n/a, n/a, n/a]
 *
 * where T is the total time of flight
 * 
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE mga_target_event: public base
{
	public:
	mga_target_event(const kep_toolbox::planet::planet_ptr start = kep_toolbox::planet::jpl_lp("earth").clone(),
						 const kep_toolbox::planet::planet_ptr end = kep_toolbox::planet::jpl_lp("mars").clone(),
						 const kep_toolbox::epoch t_end = kep_toolbox::epoch(0.0),
						 double T_max = 700.0,
						 bool discount_launcher = true
						 );

		mga_target_event(const mga_target_event&);
		base_ptr clone() const;
		
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		std::string human_readable_extra() const;
		
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_start;
			ar & m_end;
			ar & m_t_end;
			ar & m_T_max;
			ar & m_discount_launcher;
		}
	kep_toolbox::planet::planet_ptr m_start;
	kep_toolbox::planet::planet_ptr m_end;
	kep_toolbox::epoch m_t_end;
	double m_T_max;
	bool m_discount_launcher;
};

}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::mga_target_event)

#endif // PAGMO_PROBLEM_MGA_TARGET_EVENT_H
