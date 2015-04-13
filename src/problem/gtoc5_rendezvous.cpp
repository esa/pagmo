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
#include "gtoc5_rendezvous.h"
#include <keplerian_toolbox/astro_constants.h>
#include <keplerian_toolbox/planet/gtoc5.h>

using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {
/// Constructor.

gtoc5_rendezvous::gtoc5_rendezvous(int segments, int source, int target, const double &lb_epoch, const double  &initial_mass, const double &ctol):
	base(segments * 3 + 3, 0, 1, 7 + segments, segments,ctol),
	m_n_segments(segments),m_source(source),m_target(target),m_lb_epoch(lb_epoch),m_initial_mass(initial_mass)
{
	std::vector<double> lb_v(get_dimension());
	std::vector<double> ub_v(get_dimension());

	// Start (MJD).
	lb_v[0] = m_lb_epoch;
	ub_v[0] = m_lb_epoch + 356.25 * 10;

	// Leg duration.
	lb_v[1] = 10;
	ub_v[1] = 365.25 * 4;

	// Final mass.
	lb_v[2] = 500;
	ub_v[2] = m_initial_mass;

	// I Throttles
	for (int i = 3; i < segments * 3 + 3; ++i)
	{
		lb_v[i] = -1;
		ub_v[i] = 1;
	}

	//Copying the lb,ub vector into the problems bounds
	set_bounds(lb_v,ub_v);
	
	// Set spacecraft.
	m_leg.set_spacecraft(kep_toolbox::sims_flanagan::spacecraft(m_initial_mass,0.3,3000));
}

/// Clone method.
base_ptr gtoc5_rendezvous::clone() const
{
	return base_ptr(new gtoc5_rendezvous(*this));
}

/// Implementation of the objective function.
void gtoc5_rendezvous::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = -x[2];
}

/// Implementation of the constraint function.
void gtoc5_rendezvous::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_i(x[0],epoch::MJD), epoch_f(x[1] + x[0],epoch::MJD);
	array3D v0, r0, vf, rf;
	m_source.eph(epoch_i,r0,v0);
	m_target.eph(epoch_f,rf,vf);
	m_leg.set_leg(epoch_i,sc_state(r0,v0,m_leg.get_spacecraft().get_mass()),x.begin() + 3, x.end(),epoch_f,sc_state(rf,vf,x[2]),ASTRO_MU_SUN);

	// We evaluate the state mismatch at the mid-point. And we use astronomical units to scale them
	m_leg.get_mismatch_con(c.begin(), c.begin() + 7);
	for (int i=0; i<3; ++i) c[i]/=ASTRO_AU;
	for (int i=3; i<6; ++i) c[i]/=ASTRO_EARTH_VELOCITY;
	c[6] /= m_leg.get_spacecraft().get_mass();
	// We evaluate the constraints on the throttles writing on the 7th mismatch constrant (mass is off)
	m_leg.get_throttles_con(c.begin() + 7, c.begin() + 7 + m_n_segments);
}

/// Implementation of the sparsity structure: automated detection
void gtoc5_rendezvous::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
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

std::string gtoc5_rendezvous::get_name() const
{
	return "GTOC5 Rendezvous phase";
}

std::string gtoc5_rendezvous::pretty(const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_i(x[0],epoch::MJD), epoch_f(x[1] + x[0],epoch::MJD);
	array3D v0, r0, vf, rf;
	m_source.eph(epoch_i,r0,v0);
	m_target.eph(epoch_f,rf,vf);
	m_leg.set_leg(epoch_i,sc_state(r0,v0,m_leg.get_spacecraft().get_mass()),x.begin() + 3, x.end(),epoch_f,sc_state(rf,vf,x[2]),ASTRO_MU_SUN);

	std::ostringstream oss;
	oss << m_leg << '\n' << m_source << '\n' << m_target << '\n';
	return oss.str();
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::gtoc5_rendezvous);
