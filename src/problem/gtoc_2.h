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

#ifndef PAGMO_PROBLEM_GTOC_2_H
#define PAGMO_PROBLEM_GTOC_2_H

#include <vector>
#include <string>
#include <keplerian_toolbox/epoch.h>
#include <keplerian_toolbox/planet/gtoc2.h>
#include <keplerian_toolbox/sims_flanagan/leg.h>
#include <keplerian_toolbox/sims_flanagan/spacecraft.h>
#include <keplerian_toolbox/sims_flanagan/throttle.h>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"


namespace pagmo { namespace problem {


/// GTOC_2 Low-Thrust Multiple Asteroid Randezvous Problem
/**
 * This is the problem given by Jet Propulsion Laboratories as the 2nd Global Trajectory
 * Optimization Competition. It is transcribed assembling 4 Sims-Flanagan trajectories legs
 *
 * @see http://www.esa.int/gsp/ACT/inf/op/globopt/evevejsa.htm
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE gtoc_2: public base
{
	public:
		/// The objective function can be defined as final mass, final time or mass/time
		enum objective {MASS,TIME,MASS_TIME};
		/// Constructor
		gtoc_2(int = 815, int = 300, int = 110, int = 47, int = 10, objective = MASS_TIME);
		base_ptr clone() const;
		//void set_sparsity(int &, std::vector<int> &, std::vector<int> &) const;
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x) const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		std::string human_readable_extra() const;
	private:
		template <class Iterator>
		kep_toolbox::sims_flanagan::throttle get_nth_throttle(int n, Iterator it, const kep_toolbox::epoch &start, const kep_toolbox::epoch &end) const
		{
			Iterator n_it = it + 3 * n;
			kep_toolbox::array3D tmp = {{ *n_it, *(n_it + 1), *(n_it + 2)}};
			const double seg_duration = (end.mjd() - start.mjd()) / m_n_seg;
			return kep_toolbox::sims_flanagan::throttle( kep_toolbox::epoch(start.mjd() + seg_duration * n,kep_toolbox::epoch::MJD),
				kep_toolbox::epoch(start.mjd() + seg_duration * (n + 1),kep_toolbox::epoch::MJD),
				tmp);
		}
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<int &>(m_n_seg);
			ar & m_asteroids;
			ar & m_legs;
			ar & const_cast< kep_toolbox::sims_flanagan::spacecraft &>(m_spacecraft);
			ar & m_obj;
		}
		const int												m_n_seg;
		std::vector<kep_toolbox::planet::gtoc2>				m_asteroids;
		mutable std::vector< kep_toolbox::sims_flanagan::leg>	m_legs;
		const kep_toolbox::sims_flanagan::spacecraft			m_spacecraft;
		objective												m_obj;
};

} } // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::gtoc_2)

#endif

