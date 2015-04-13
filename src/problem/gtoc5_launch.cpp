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
#include "gtoc5_launch.h"
#include <keplerian_toolbox/sims_flanagan/codings.h>
#include <keplerian_toolbox/astro_constants.h>
#include <keplerian_toolbox/planet/gtoc5.h>

using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {
/// Constructor.

gtoc5_launch::gtoc5_launch(int segments, int target, objective obj, const double &ctol) :
	base(segments * 3 + 6, 0, 1, 7 + segments + 1, segments + 1,ctol),
	m_n_segments(segments),m_earth(),m_target(target),m_obj(obj)
{
	std::vector<double> lb_v(get_dimension());
	std::vector<double> ub_v(get_dimension());

	// Start (MJD).
	lb_v[0] = 57023;
	ub_v[0] = 61041;

	// Leg duration in days
	lb_v[1] = 30;
	ub_v[1] = 365.25 * 3;

	// Final mass.
	lb_v[2] = 500;
	ub_v[2] = 4000;

	// Start Velocity
	lb_v[3] = -5000;
	lb_v[4] = -5000;
	lb_v[5] = -5000;
	ub_v[3] = 5000;
	ub_v[4] = 5000;
	ub_v[5] = 5000;

	// I Throttles
	for (int i = 6; i < segments * 3 + 6; ++i)
	{
		lb_v[i] = -1;
		ub_v[i] = 1;
	}

	//Copying the lb,ub vector into the problems bounds
	set_bounds(lb_v,ub_v);
	
	// Set GTOC5 spacecraft.
	m_leg.set_spacecraft(kep_toolbox::sims_flanagan::spacecraft(4000,0.3,3000));
}

/// Clone method.
base_ptr gtoc5_launch::clone() const
{
	return base_ptr(new gtoc5_launch(*this));
}

/// Implementation of the objective function.
void gtoc5_launch::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	if (m_obj == MASS) {
		f[0] = -x[2] / m_leg.get_spacecraft().get_mass();
	} else {
		f[0] = x[1] / 365.25;
	}
}

/// Implementation of the constraint function.
void gtoc5_launch::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	using namespace kep_toolbox;
	// 1 - We set the leg.
	const epoch epoch_i(x[0],epoch::MJD), epoch_f(x[1] + x[0],epoch::MJD);
	array3D v0, r0, vf, rf;
	m_earth.eph(epoch_i,r0,v0);
	m_target.eph(epoch_f,rf,vf);

	v0[0] += x[3];
	v0[1] += x[4];
	v0[2] += x[5];
	m_leg.set_leg(epoch_i,sc_state(r0,v0,m_leg.get_spacecraft().get_mass()),x.begin() + 6, x.end(),epoch_f,sc_state(rf,vf,x[2]),ASTRO_MU_SUN);

	// We evaluate the state mismatch at the mid-point. And we use astronomical units to scale them
	m_leg.get_mismatch_con(c.begin(), c.begin() + 7);
	for (int i=0; i<3; ++i) c[i]/=ASTRO_AU;
	for (int i=3; i<6; ++i) c[i]/=ASTRO_EARTH_VELOCITY;
	c[6] /= m_leg.get_spacecraft().get_mass();
	// We evaluate the constraints on the throttles writing on the 7th mismatch constrant (mass is off)
	m_leg.get_throttles_con(c.begin() + 7, c.begin() + 7 + m_n_segments);
	c[7 + m_n_segments] = (x[3]*x[3] + x[4]*x[4] + x[5]*x[5] - 25000000) / ASTRO_EARTH_VELOCITY / ASTRO_EARTH_VELOCITY;
}

/// Implementation of the sparsity structure: automated detection
void gtoc5_launch::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
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

std::string gtoc5_launch::get_name() const
{
	return "GTOC5 Launch phase";
}

std::string gtoc5_launch::pretty(const decision_vector &x) const
{
	using namespace kep_toolbox;
	// 1 - We set the leg.
	const epoch epoch_i(x[0],epoch::MJD), epoch_f(x[1] + x[0],epoch::MJD);
	array3D v0, r0, vf, rf;
	m_earth.eph(epoch_i,r0,v0);
	m_target.eph(epoch_f,rf,vf);

	v0[0] += x[2];
	v0[1] += x[3];
	v0[2] += x[4];
	m_leg.set_leg(epoch_i,sc_state(r0,v0,m_leg.get_spacecraft().get_mass()),x.begin() + 6, x.end(),epoch_f,sc_state(rf,vf,x[5]),ASTRO_MU_SUN);

	std::ostringstream oss;
	oss << m_leg << '\n' << m_earth << '\n' << m_target << '\n';
	return oss.str();
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::gtoc5_launch);
