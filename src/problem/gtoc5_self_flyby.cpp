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
#include "gtoc5_self_flyby.h"
#include <keplerian_toolbox/astro_constants.h>
#include <keplerian_toolbox/planet/gtoc5.h>

using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {
/// Constructor.

gtoc5_self_flyby::gtoc5_self_flyby(int segments, int ast_id, const double &mjd, const double  &initial_mass, const double &ctol):
	base(segments * 3 + 5, 0, 1, 7 + segments + 1, segments + 1,ctol),
	m_n_segments(segments),m_ast(ast_id),m_mjd(mjd),m_initial_mass(initial_mass)
{
	std::vector<double> lb_v(get_dimension());
	std::vector<double> ub_v(get_dimension());

	// Flyby (MJD) transfer time in days
	lb_v[0] = 5;
	ub_v[0] = 200;

	// Mass at fly-by
	lb_v[1] = 100;
	ub_v[1] = m_initial_mass;

	// Velocity at fly-by
	for (int i = 2; i < 5; ++i) {
		lb_v[i] = -1500;
		ub_v[i] = 1500;
	}

	// I Throttles
	for (int i = 5; i < segments * 3 + 5; ++i)
	{
		lb_v[i] = -1;
		ub_v[i] = 1;
	}

	set_bounds(lb_v,ub_v);
	
	// Set spacecraft.
	m_leg.set_spacecraft(kep_toolbox::sims_flanagan::spacecraft(m_initial_mass,0.3,3000));
}

/// Clone method.
base_ptr gtoc5_self_flyby::clone() const
{
	return base_ptr(new gtoc5_self_flyby(*this));
}

/// Implementation of the objective function.
void gtoc5_self_flyby::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = x[0];   //Time of Flight
}

/// Implementation of the constraint function.
void gtoc5_self_flyby::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_source(m_mjd,epoch::MJD), epoch_target(m_mjd + x[0],epoch::MJD);
	array3D v_source, r_source, v_target, r_target;
	m_ast.eph(epoch_source,r_source,v_source);
	m_ast.eph(epoch_target,r_target,v_target);
	v_target[0] += x[2];
	v_target[1] += x[3];
	v_target[2] += x[4];
	m_leg.set_leg(epoch_source,sc_state(r_source,v_source,m_leg.get_spacecraft().get_mass()),x.begin() + 5, x.begin() + 5 + m_n_segments * 3,
		epoch_target,sc_state(r_target,v_target,x[1]),ASTRO_MU_SUN);


	// We evaluate the state mismatch at the mid-point. And we use astronomical units to scale them
	m_leg.get_mismatch_con(c.begin(), c.begin() + 7);
	for (int i=0; i<3; ++i) c[i]/=ASTRO_AU;
	for (int i=3; i<6; ++i) c[i]/=ASTRO_EARTH_VELOCITY;
	c[6] /= m_leg.get_spacecraft().get_mass();

	// We evaluate the constraints on the throttles writing on the 7th mismatch constrant (mass is off)
	m_leg.get_throttles_con(c.begin() + 7, c.begin() + 7 + m_n_segments);
	
	// Minimum flyby speed = 0.4 km/s.
	c[7 + m_n_segments] = (160000. - (x[2] * x[2] + x[3] * x[3] + x[4] * x[4]));
}

/// Implementation of the sparsity structure: automated detection
void gtoc5_self_flyby::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
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

std::string gtoc5_self_flyby::get_name() const
{
	return "GTOC5 Flyby phase";
}

std::string gtoc5_self_flyby::pretty(const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_source(m_mjd,epoch::MJD), epoch_target(m_mjd + x[0],epoch::MJD);
	array3D v_source, r_source, v_target, r_target;
	m_ast.eph(epoch_source,r_source,v_source);
	m_ast.eph(epoch_target,r_target,v_target);
	v_target[0] += x[2];
	v_target[1] += x[3];
	v_target[2] += x[4];
	m_leg.set_leg(epoch_source,sc_state(r_source,v_source,m_leg.get_spacecraft().get_mass()),x.begin() + 5, x.begin() + 5 + m_n_segments * 3,
		epoch_target,sc_state(r_target,v_target,x[1]),ASTRO_MU_SUN);

	std::ostringstream oss;
	oss << m_leg << '\n' << m_ast << '\n';
	return oss.str();
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::gtoc5_self_flyby);
