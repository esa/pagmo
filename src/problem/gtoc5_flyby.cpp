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
#include "gtoc5_flyby.h"
#include <keplerian_toolbox/astro_constants.h>
#include <keplerian_toolbox/planet/gtoc5.h>

using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {
/// Constructor.

gtoc5_flyby::gtoc5_flyby(int segments, int source, int flyby, int target, const double &lb_epoch, const double  &initial_mass, objective obj, const double & tof_ub, const double &ctol):
	base(segments * 6 + 8, 0, 1, 14 + 2 * segments + 1, 2 * segments + 1,ctol),
	m_n_segments(segments),m_source(source),m_flyby(flyby),m_target(target),m_lb_epoch(lb_epoch),m_initial_mass(initial_mass),m_obj(obj)
{
	std::vector<double> lb_v(get_dimension());
	std::vector<double> ub_v(get_dimension());

	// Source (MJD).
	lb_v[0] = m_lb_epoch;
	ub_v[0] = m_lb_epoch + 200;

	// Flyby days fraction
	lb_v[1] = 0.001;
	ub_v[1] = 0.99;

	// Total transfer time in days
	lb_v[2] = 10;
	ub_v[2] = 365.25 * tof_ub;

	// Mass at fly-by
	lb_v[3] = 5;
	ub_v[3] = m_initial_mass;

	// Mass at target
	lb_v[4] = 5;
	ub_v[4] = m_initial_mass;

	// Velocity at fly-by
	for (int i = 5; i < 8; ++i) {
		lb_v[i] = -2000;
		ub_v[i] = 2000;
	}

	// I Throttles
	for (int i = 8; i < segments * 6 + 8; ++i)
	{
		lb_v[i] = -1;
		ub_v[i] = 1;
	}

	set_bounds(lb_v,ub_v);
	
	// Set spacecraft.
	m_leg1.set_spacecraft(kep_toolbox::sims_flanagan::spacecraft(m_initial_mass,0.3,3000));
	m_leg2.set_spacecraft(kep_toolbox::sims_flanagan::spacecraft(m_initial_mass,0.3,3000));
}

/// Clone method.
base_ptr gtoc5_flyby::clone() const
{
	return base_ptr(new gtoc5_flyby(*this));
}

/// Implementation of the objective function.
void gtoc5_flyby::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	if (m_obj == MASS) {
		f[0] = -x[4];
	} else if (m_obj == TIME) {
		f[0] = x[2];
	} else {
		f[0] = x[0] + x[2];
	}
}

/// Implementation of the constraint function.
void gtoc5_flyby::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_source(x[0],epoch::MJD), epoch_flyby(x[0] + x[1]*x[2],epoch::MJD), epoch_target(x[0] + x[2],epoch::MJD);
	array3D v_source, r_source, v_flyby, r_flyby, v_target, r_target;
	m_source.eph(epoch_source,r_source,v_source);
	m_flyby.eph(epoch_flyby,r_flyby,v_flyby);
	v_flyby[0] += x[5];
	v_flyby[1] += x[6];
	v_flyby[2] += x[7];
	m_target.eph(epoch_target,r_target,v_target);
	m_leg1.set_leg(epoch_source,sc_state(r_source,v_source,m_leg1.get_spacecraft().get_mass()),x.begin() + 8, x.begin() + 8 + m_n_segments * 3,
		epoch_flyby,sc_state(r_flyby,v_flyby,x[3]),ASTRO_MU_SUN);
	m_leg2.set_leg(epoch_flyby,sc_state(r_flyby,v_flyby,x[3] - 1),x.begin() + 8 + m_n_segments * 3, x.end(),
		epoch_target,sc_state(r_target,v_target,x[4]),ASTRO_MU_SUN);

	// We evaluate the state mismatch at the mid-point. And we use astronomical units to scale them
	m_leg1.get_mismatch_con(c.begin(), c.begin() + 7);
	for (int i=0; i<3; ++i) c[i]/=ASTRO_AU;
	for (int i=3; i<6; ++i) c[i]/=ASTRO_EARTH_VELOCITY;
	c[6] /= m_leg1.get_spacecraft().get_mass();

	m_leg2.get_mismatch_con(c.begin() + 7, c.begin() + 14);
	for (int i=7; i<10; ++i) c[i]/=ASTRO_AU;
	for (int i=10; i<13; ++i) c[i]/=ASTRO_EARTH_VELOCITY;
	c[13] /= m_leg1.get_spacecraft().get_mass(); // Scaling units for mass are identical throughtout the legs, thus m_leg1.

	// We evaluate the constraints on the throttles writing on the 7th mismatch constrant (mass is off)
	m_leg1.get_throttles_con(c.begin() + 14, c.begin() + 14 + m_n_segments);
	m_leg2.get_throttles_con(c.begin() + 14 + m_n_segments, c.begin() + 14 + 2 * m_n_segments);
	
	// Minimum flyby speed = 0.4 km/s.
	c[14 + 2 * m_n_segments] = (160000. - (x[5] * x[5] + x[6] * x[6] + x[7] * x[7])) / 160000.;
}

/// Implementation of the sparsity structure: automated detection
void gtoc5_flyby::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
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

std::string gtoc5_flyby::get_name() const
{
	return "GTOC5 Flyby phase";
}

std::string gtoc5_flyby::pretty(const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_source(x[0],epoch::MJD), epoch_flyby(x[0] + x[1],epoch::MJD), epoch_target(x[0] + x[1] + x[2],epoch::MJD);
	array3D v_source, r_source, v_flyby, r_flyby, v_target, r_target;
	m_source.eph(epoch_source,r_source,v_source);
	m_flyby.eph(epoch_flyby,r_flyby,v_flyby);
	v_flyby[0] += x[5];
	v_flyby[1] += x[6];
	v_flyby[2] += x[7];
	m_target.eph(epoch_target,r_target,v_target);
	m_leg1.set_leg(epoch_source,sc_state(r_source,v_source,m_leg1.get_spacecraft().get_mass()),x.begin() + 8, x.begin() + 8 + m_n_segments * 3,
		epoch_flyby,sc_state(r_flyby,v_flyby,x[3]),ASTRO_MU_SUN);
	m_leg2.set_leg(epoch_flyby,sc_state(r_flyby,v_flyby,x[3] - 1),x.begin() + 8 + m_n_segments * 3, x.end(),
		epoch_target,sc_state(r_target,v_target,x[4]),ASTRO_MU_SUN);

	std::ostringstream oss;
	oss << m_leg1 << '\n' << m_source << '\n' << m_flyby << '\n';
	oss << m_leg2 << '\n' << m_flyby << '\n' << m_target << '\n';
	return oss.str();
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::gtoc5_flyby);
