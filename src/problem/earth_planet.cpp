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

#include <boost/integer_traits.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "earth_planet.h"
#include "../keplerian_toolbox/codings.h"
#include "../keplerian_toolbox/planet.h"

using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {
/// Constructor.

earth_planet::earth_planet(int segments, planet::common_name planet_id) : base(base_format(1,segments,1000).size(), 0, 1, 6 + segments + 1 +1, segments+2),
 encoding(1,segments,1000), vmax(3000),n_segments(segments)
{
	std::vector<double> lb_v(get_dimension());
	std::vector<double> ub_v(get_dimension());

	//Start
	lb_v[encoding.leg_start_epoch_i(0)[0]] = 0;
	ub_v[encoding.leg_start_epoch_i(0)[0]] = 1000;

	//End
	lb_v[encoding.leg_end_epoch_i(0)[0]] = 500;
	ub_v[encoding.leg_end_epoch_i(0)[0]] = 1500;

	//Start Velocity
	std::vector<int> tmp = encoding.leg_start_velocity_i(0);
	for (std::vector<int>::size_type i = 0; i < tmp.size() ; ++i)
	{
		lb_v[tmp[i]] = -3;
		ub_v[tmp[i]] = 3;
	}

	//End Velocity
	tmp = encoding.leg_end_velocity_i(0);
	for (std::vector<int>::size_type i = 0; i < tmp.size() ; ++i)
	{
		lb_v[tmp[i]] = 0;
		ub_v[tmp[i]] = 0;
	}

	//I Throttles
	for (int j = 0; j<encoding.n_segments(0); ++j)
	{
		tmp = encoding.segment_thrust_i(0,j);
		for (std::vector<int>::size_type i = 0; i < tmp.size() ; ++i)
		{
			lb_v[tmp[i]] = 0;
			ub_v[tmp[i]] = 1;
		}
	}

	//traj_fb constructor
	std::vector<planet> sequence;
	sequence.push_back(planet(planet::EARTH));
	sequence.push_back(planet(planet_id));
	trajectory = fb_traj(sequence,segments,1000,0.05,4500);
}

/// Clone method.
base_ptr earth_planet::clone() const
{
	return base_ptr(new earth_planet(*this));
}

/// Implementation of the objective function.
void earth_planet::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	trajectory.init_from_full_vector(x.begin(),x.end(),encoding);
	f[0] = trajectory.get_leg(0).evaluate_dv();
}

/// Implementation of the constraint function.
void earth_planet::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	trajectory.init_from_full_vector(x.begin(),x.end(),encoding);
	trajectory.evaluate_all_mismatch_con(c.begin(), c.begin() + 7);
	trajectory.get_leg(0).get_throttles_con(c.begin() + 6, c.begin() + 6 + n_segments);\

	c[6 + n_segments] = trajectory.evaluate_leg_vinf2_i(0) - vmax*vmax;
	c[7 + n_segments] = trajectory.get_leg(0).get_t_i().mjd2000() - trajectory.get_leg(0).get_t_f().mjd2000();
}

/// Implementation of the sparsity structure: automated detection
void earth_planet::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	//Initial point
	decision_vector x0(get_dimension(),1);
	//Numerical procedure
	estimate_sparsity(x0, lenG, iGfun, jGvar);
}

} }
