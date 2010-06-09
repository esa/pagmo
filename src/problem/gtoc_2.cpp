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

#include <algorithm>
#include <functional>
#include <vector>
#include <string>
#include <sstream>

#include "../exceptions.h"
#include "../keplerian_toolbox/keplerian_toolbox.h"
#include "../types.h"
#include "base.h"
#include "gtoc_2.h"

using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {

gtoc_2::gtoc_2(int ast1, int ast2, int ast3, int ast4, int n_seg):
	base(12 * n_seg + 15,0,1,7 * 4 + n_seg * 4 + 1, n_seg * 4 + 1, 1E-9),
	m_n_seg(n_seg),m_spacecraft(1500.,.1,4000.)
{
	if (n_seg <= 0) {
		pagmo_throw(value_error,"invalid number of segments");
	}
	m_asteroids.push_back(asteroid_gtoc2(910));
	m_asteroids.push_back(asteroid_gtoc2(ast1));
	m_asteroids.push_back(asteroid_gtoc2(ast2));
	m_asteroids.push_back(asteroid_gtoc2(ast3));
	m_asteroids.push_back(asteroid_gtoc2(ast4));
	// Build legs.
	for (int i = 0; i < 4; ++i) {
		m_legs.push_back(leg());
		m_legs.back().set_spacecraft(m_spacecraft);
		m_legs.back().set_mu(1.32712440018e20);
		m_legs.back().set_throttles_size(n_seg);
	}
	// Check that input asteroids belong to different groups.
	int asteroid_groups[4] = {m_asteroids[1].get_group(), m_asteroids[2].get_group(),
		m_asteroids[3].get_group(), m_asteroids[4].get_group()};
	std::sort(asteroid_groups,asteroid_groups + 4);
	if (std::unique(asteroid_groups,asteroid_groups + 4) - asteroid_groups != 4) {
		pagmo_throw(value_error,"asteroids must belong to different groups");
	}
	if (std::find(asteroid_groups,asteroid_groups + 4,5) != asteroid_groups + 4) {
		pagmo_throw(value_error,"the Earth cannot appear in the list of asteroids");
	}
	decision_vector lb_v, ub_v;
	// Launch date.
	lb_v.push_back(57023.5);
	ub_v.push_back(64693.5);
	// Flight and waiting times.
	for (int i = 0; i < 3; ++i) {
		// Flight time.
		lb_v.push_back(5);
		ub_v.push_back(1800);
		// Waiting time.
		lb_v.push_back(90);
		ub_v.push_back(500);
	}
	// Final transfer time.
	lb_v.push_back(5);
	ub_v.push_back(1000);
	// Masses.
	for (int i = 0; i < 4; ++i) {
		lb_v.push_back(500);
		ub_v.push_back(1500);
	}
	// Vinf cartesian km/s.
	for (int i = 0; i < 3; ++i) {
		lb_v.push_back(-3.5);
		ub_v.push_back(3.5);
	}
	// Throttles.
	for (int i = 0; i < 12 * n_seg; ++i) {
		lb_v.push_back(-1);
		ub_v.push_back(1);
	}
	set_bounds(lb_v,ub_v);
}

base_ptr gtoc_2::clone() const
{
	return base_ptr(new gtoc_2(*this));
}

void gtoc_2::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	// Objective function is m_f / t_f.
	f[0] = x[11] / (std::accumulate(x.begin() + 1,x.begin() + 8,0) * ASTRO_DAY2YEAR);
}



void gtoc_2::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	// Cached values.
	array3D r, v;
	double initial_mass = m_spacecraft.get_mass();
	double final_mass = x[8];
	// Build legs.
	for (int i = 0; i < 4; ++i) {
		// Calculate start-end leg epochs.
		epoch start(std::accumulate(x.begin(),x.begin() + 2 * i + 1, 0),epoch::MJD),
			end(std::accumulate(x.begin(),x.begin() + 2 * i + 2, 0),epoch::MJD);
		// Set leg's start-end epochs.
		m_legs[i].set_t_i(start);
		m_legs[i].set_t_f(end);
		// Initial state.
		m_asteroids[i].get_eph(start,r,v);
		// First leg: Vinf is added to the state.
		if (i == 0) {
			for (int j = 0; j < 3; ++j) {
				v[j] += x[12 + j] * 1000;
			}
		}
		m_legs[i].set_x_i(sc_state(r,v,initial_mass));
		// Throttles.
		for (int j = 0; j < m_n_seg; ++j) {
			m_legs[i].set_throttle(j,get_nth_throttle(j,x.begin() + 15 + 3 * i * m_n_seg,start,end));
		}
		// Final state.
		m_asteroids[i + 1].get_eph(end,r,v);
		m_legs[i].set_x_f(sc_state(r,v,final_mass));
		// Update masses.
		initial_mass = final_mass;
		// Do not update mass if final iteration (useless).
		if (i != 3) { 
			final_mass = x[9 + i];
		}
	}
	// Load state mismatches into constraints vector.
	for (int i = 0; i < 4; ++i) {
		m_legs[i].get_mismatch_con(c.begin() + 7 * i, c.begin() + 7 * (i + 1));
		// Passing non-dimensional units to the solver.
		for (int j = 0; j < 3; ++j) {
			c[7 * i + j] /= ASTRO_AU;
			c[7 * i + j + 3] /= ASTRO_EARTH_VELOCITY;
		}
	}
	// Throttles constraints.
	for (int i = 0; i < 4; ++i) {
		m_legs[i].get_throttles_con(c.begin() + 28 + i * m_n_seg, c.begin() + 28 + (i + 1) * m_n_seg);
	}
	// Vinf constraint.
	c.back() = (x[12] * x[12] + x[13] * x[13] + x[14] * x[14] - 3.5 * 3.5) / (ASTRO_EARTH_VELOCITY * ASTRO_EARTH_VELOCITY) * 1000000;
}

/// Implementation of the sparsity structure: automated detection
// void gtoc_2::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
// {
// 	//Initial point
// 	decision_vector x0(get_dimension());
// 	for (pagmo::decision_vector::size_type i = 0; i<x0.size(); ++i)
// 	{
// 		x0[i] = get_lb()[i] + (get_ub()[i] - get_lb()[i]) / 3.12345;
// 	}
// 	//Numerical procedure
// 	estimate_sparsity(x0, lenG, iGfun, jGvar);
// }

std::string gtoc_2::get_name() const
{
	return "GTOC_2";
}

std::string gtoc_2::human_readable_extra() const {
	std::ostringstream oss;
	oss << "Asteroid Sequence: ";
	for (int i =1 ; i< 5; ++i ) oss << m_asteroids[i].get_name() << " ";
	oss << std::endl;
	return oss.str();
}

}} // namespaces
