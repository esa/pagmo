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

#include <string>

#include <boost/integer_traits.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "earth_gtoc5_asteroid.h"
#include "../keplerian_toolbox/sims_flanagan/codings.h"
#include "../keplerian_toolbox/astro_constants.h"
#include "../keplerian_toolbox/asteroid_gtoc5.h"

using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {
/// Constructor.

earth_gtoc5_asteroid::earth_gtoc5_asteroid(int segments, int target, const double &ctol) : base(segments * 3 + 5, 0, 1, 7 + segments + 1, segments + 1,ctol),
 m_n_segments(segments),m_earth(),m_target(target)
{
	std::vector<double> lb_v(get_dimension());
	std::vector<double> ub_v(get_dimension());

	// Start (MJD).
	lb_v[0] = 57023;
	ub_v[0] = 61041;

	// End (MJD)
	lb_v.back() = 57023;
	ub_v.back() = 61041 + 15 * 365.25;

	// Start Velocity orientation (magnitude fixed to 5 km/s) -> sphere picking.
	lb_v[1] = -10;
	lb_v[2] = -10;
	ub_v[1] = 10;
	ub_v[2] = 10;

	// I Throttles
	for (int i = 3; i < segments * 3 + 3; ++i)
	{
		lb_v[i] = -1;
		ub_v[i] = 1;
	}

	// Final mass.
	lb_v[get_dimension() - 2] = 500;
	ub_v[get_dimension() - 2] = 4000;

	set_bounds(lb_v,ub_v);
	
	// Set spacecraft.
	m_leg.set_spacecraft(kep_toolbox::spacecraft(4000,0.3,3000));
}

/// Clone method.
base_ptr earth_gtoc5_asteroid::clone() const
{
	return base_ptr(new earth_gtoc5_asteroid(*this));
}

/// Implementation of the objective function.
void earth_gtoc5_asteroid::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = -x[get_dimension() - 2];
}

/// Implementation of the constraint function.
void earth_gtoc5_asteroid::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_i(x[0],epoch::MJD), epoch_f(x.back(),epoch::MJD);
	array3D v0, r0, vf, rf;
	m_earth.get_eph(epoch_i,r0,v0);
	m_target.get_eph(epoch_f,rf,vf);
	const double theta = x[1], phi = x[2];
	v0[0] += 5000 * std::cos(theta) * std::sin(phi);
	v0[1] += 5000 * std::sin(theta) * std::sin(phi);
	v0[2] += 5000 * std::cos(phi);
	m_leg.set_leg(epoch_i,sc_state(r0,v0,m_leg.get_spacecraft().get_mass()),x.begin() + 3, x.end() - 2,epoch_f,sc_state(rf,vf,x[get_dimension() - 2]),ASTRO_MU_SUN);

	// We evaluate the state mismatch at the mid-point. And we use astronomical units to scale them
	m_leg.get_mismatch_con(c.begin(), c.begin() + 7);
	for (int i=0; i<3; ++i) c[i]/=ASTRO_AU;
	for (int i=3; i<6; ++i) c[i]/=ASTRO_EARTH_VELOCITY;
	// We evaluate the constraints on the throttles writing on the 7th mismatch constrant (mass is off)
	m_leg.get_throttles_con(c.begin() + 7, c.begin() + 7 + m_n_segments);
	// We evaluate the linear constraint on the epochs (tf > ti)
	c[7 + m_n_segments] = x.back() - x.front();
}

/// Implementation of the sparsity structure: automated detection
void earth_gtoc5_asteroid::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	//Initial point
	decision_vector x0(get_dimension());
	for (pagmo::decision_vector::size_type i = 0; i<x0.size(); ++i)
	{
		x0[i] = get_lb()[i] + (get_ub()[i] - get_lb()[i]) / 3.12345;
	}
	//Numerical procedure
	estimate_sparsity(x0, lenG, iGfun, jGvar);
}

std::string earth_gtoc5_asteroid::get_name() const
{
	return "Earth-GTOC5-Asteroid";
}

std::string earth_gtoc5_asteroid::pretty(const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_i(x[0],epoch::MJD), epoch_f(x.back(),epoch::MJD);
	array3D v0, r0, vf, rf;
	m_earth.get_eph(epoch_i,r0,v0);
	m_target.get_eph(epoch_f,rf,vf);
	const double theta = x[1], phi = x[2];
	v0[0] += 5000 * std::cos(theta) * std::sin(phi);
	v0[1] += 5000 * std::sin(theta) * std::sin(phi);
	v0[2] += 5000 * std::cos(phi);
	m_leg.set_leg(epoch_i,sc_state(r0,v0,m_leg.get_spacecraft().get_mass()),x.begin() + 3, x.end() - 2,epoch_f,sc_state(rf,vf,x[get_dimension() - 2]),ASTRO_MU_SUN);
	
	std::ostringstream oss;
	oss << m_leg << '\n' << m_earth << '\n' << m_target << '\n';
	return oss.str();
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::earth_gtoc5_asteroid);
