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

#ifndef PAGMO_GTOC5_RENDEZVOUS_H
#define PAGMO_GTOC5_RENDEZVOUS_H

#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include <keplerian_toolbox/keplerian_toolbox.h>
#include "base.h"

namespace pagmo { namespace problem {

/// Test problem kep tool
/**
 *
 *
 * @author Dario Izzo (dario.izzo@esa.int)
 */

class __PAGMO_VISIBLE gtoc5_rendezvous: public base
{
	public:
		gtoc5_rendezvous(int = 10, int = 1, int = 2, const double & = 57023, const double & = 4000, const double & = 1E-5);
		base_ptr clone() const;
		std::string get_name() const;
		/// A nice string representation of a chromosome
		std::string pretty(const decision_vector &) const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		void set_sparsity(int &, std::vector<int> &, std::vector<int> &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_n_segments;
			ar & const_cast<kep_toolbox::planet::gtoc5 &>(m_source);
			ar & const_cast<kep_toolbox::planet::gtoc5 &>(m_target);
			ar & const_cast<double &>(m_lb_epoch);
			ar & const_cast<double &>(m_initial_mass);
			ar & m_leg;
		}
		int 						m_n_segments;
		const kep_toolbox::planet::gtoc5 		m_source;
		const kep_toolbox::planet::gtoc5 		m_target;
		const double					m_lb_epoch;
		const double					m_initial_mass;
		mutable kep_toolbox::sims_flanagan::leg		m_leg;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::gtoc5_rendezvous);

#endif // GTOC5_GTOC5_ASTEROID_H
