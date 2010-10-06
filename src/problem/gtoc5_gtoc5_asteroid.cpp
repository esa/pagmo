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
#include "gtoc5_gtoc5_asteroid.h"
#include "../keplerian_toolbox/astro_constants.h"
#include "../keplerian_toolbox/asteroid_gtoc5.h"

using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {
/// Constructor.

gtoc5_gtoc5_asteroid::gtoc5_gtoc5_asteroid(int segments, int source, int target, const double &lb_epoch, const double  &initial_mass, const double &ctol):
	base(segments * 3 + 3, 0, 1, 7 + segments + 2, segments + 2,ctol),
	m_n_segments(segments),m_source(source),m_target(target),m_lb_epoch(lb_epoch),m_initial_mass(initial_mass)
{
	std::vector<double> lb_v(get_dimension());
	std::vector<double> ub_v(get_dimension());

	// Start (MJD).
	lb_v[0] = m_lb_epoch;
	ub_v[0] = m_lb_epoch + 356.25 * 10;

	// End (MJD)
	lb_v.back() = m_lb_epoch;
	ub_v.back() = m_lb_epoch + 356.25 * 15;

	// I Throttles
	for (int i = 1; i < segments * 3 + 1; ++i)
	{
		lb_v[i] = -1;
		ub_v[i] = 1;
	}

	// Final mass.
	lb_v[get_dimension() - 2] = 500;
	ub_v[get_dimension() - 2] = m_initial_mass;

	set_bounds(lb_v,ub_v);
	
	// Set spacecraft.
	m_leg.set_spacecraft(kep_toolbox::spacecraft(m_initial_mass,0.3,3000));
}

/// Clone method.
base_ptr gtoc5_gtoc5_asteroid::clone() const
{
	return base_ptr(new gtoc5_gtoc5_asteroid(*this));
}

/// Implementation of the objective function.
void gtoc5_gtoc5_asteroid::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = -x[get_dimension() - 2];
}

/// Implementation of the constraint function.
void gtoc5_gtoc5_asteroid::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_i(x[0],epoch::MJD), epoch_f(x.back(),epoch::MJD);
	array3D v0, r0, vf, rf;
	m_source.get_eph(epoch_i,r0,v0);
	m_target.get_eph(epoch_f,rf,vf);
	m_leg.set_leg(epoch_i,sc_state(r0,v0,m_leg.get_spacecraft().get_mass()),x.begin() + 1, x.end() - 2,epoch_f,sc_state(rf,vf,x[get_dimension() - 2]),ASTRO_MU_SUN);

	// We evaluate the state mismatch at the mid-point. And we use astronomical units to scale them
	m_leg.get_mismatch_con(c.begin(), c.begin() + 7);
	for (int i=0; i<3; ++i) c[i]/=ASTRO_AU;
	for (int i=3; i<6; ++i) c[i]/=ASTRO_EARTH_VELOCITY;
	c[6] /= m_leg.get_spacecraft().get_mass();
	// We evaluate the constraints on the throttles writing on the 7th mismatch constrant (mass is off)
	m_leg.get_throttles_con(c.begin() + 7, c.begin() + 7 + m_n_segments);
	// We evaluate the linear constraint on the epochs (tf > ti)
	c[7 + m_n_segments] = x.front() - x.back();
	c[8 + m_n_segments] = x.back() - x.front() - 365.25 * 5;
}

/// Implementation of the sparsity structure: automated detection
void gtoc5_gtoc5_asteroid::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
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

std::string gtoc5_gtoc5_asteroid::get_name() const
{
	return "GTOC5-GTOC5-Asteroid";
}

std::string gtoc5_gtoc5_asteroid::pretty(const decision_vector &x) const
{
	using namespace kep_toolbox;
	// We set the leg.
	const epoch epoch_i(x[0],epoch::MJD), epoch_f(x.back(),epoch::MJD);
	array3D v0, r0, vf, rf;
	m_source.get_eph(epoch_i,r0,v0);
	m_target.get_eph(epoch_f,rf,vf);
	m_leg.set_leg(epoch_i,sc_state(r0,v0,m_leg.get_spacecraft().get_mass()),x.begin() + 1, x.end() - 2,epoch_f,sc_state(rf,vf,x[get_dimension() - 2]),ASTRO_MU_SUN);
	
	std::ostringstream oss;
	oss << m_leg << '\n' << m_source << '\n' << m_target << '\n';
	return oss.str();
}

// Do just the norm.
// bool earth_gtoc5_asteroid::compare_constraints_impl(const constraint_vector &c1, const constraint_vector &c2) const
// {
// 	pagmo_assert(c1.size() == c2.size() && c1.size() == m_c_dimension);
// 	// L2 norm of constraints mismatches.
// 	double norm1 = 0, norm2 = 0;
// 	// Equality constraints.
// 	for (c_size_type i = 0; i < get_c_dimension() - get_ic_dimension(); ++i) {
// 		norm1 += std::abs(c1[i]) * std::abs(c1[i]);
// 		norm2 += std::abs(c2[i]) * std::abs(c2[i]);
// 	}
// 	// Inequality constraints.
// 	for (c_size_type i = get_c_dimension() - get_ic_dimension(); i < get_c_dimension(); ++i) {
// 		if (!test_constraint(c1,i)) {
// 			norm1 += c1[i] * c1[i];
// 		}
// 		if (!test_constraint(c2,i)) {
// 			norm2 += c2[i] * c2[i];
// 		}
// 	}
// 	return (norm1 < norm2);
// }

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::gtoc5_gtoc5_asteroid);
