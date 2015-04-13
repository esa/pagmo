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
#include <keplerian_toolbox/astro_constants.h>

#include "gtoc_2.h"
#include "../exceptions.h"
#include "../types.h"
#include "base.h"


using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {

/// Problem Constructor
/**
 * Constructs a gtoc_2 problem from an asteroid sequence
 *
 * @param[in] ast1 id of the first asteroid to visit
 * @param[in] ast2 id of the second asteroid to visit
 * @param[in] ast3 id of the third asteroid to visit
 * @param[in] ast4 id of the fourth asteroid to visit
 * @param[in] n_seg number of segments to be used per leg
 * @param[in] obj objective function in the enum {MASS,TIME,MASS_TIME}
 *
 * @throws value_error if the one or more asteroid belongs to the same group or if the
 * selected number of segments is negative
 *
 * @see problem::base constructors.
 */

gtoc_2::gtoc_2(int ast1, int ast2, int ast3, int ast4, int n_seg, objective obj):
	base(12 * n_seg + 15,0,1,7 * 4 + n_seg * 4 + 1, n_seg * 4 + 1, 1E-3),
	m_n_seg(n_seg),m_spacecraft(1500.,.1,4000.),m_obj(obj)
{
	if (n_seg <= 0) {
		pagmo_throw(value_error,"invalid number of segments");
	}
	m_asteroids.push_back(planet::gtoc2(910));
	m_asteroids.push_back(planet::gtoc2(ast1));
	m_asteroids.push_back(planet::gtoc2(ast2));
	m_asteroids.push_back(planet::gtoc2(ast3));
	m_asteroids.push_back(planet::gtoc2(ast4));
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
		ub_v.push_back(3600);
		// Waiting time.
		lb_v.push_back(90.0);
		ub_v.push_back(90.0000001);
	}
	// Final transfer time.
	lb_v.push_back(5);
	ub_v.push_back(3600);
	// Masses.
	for (int i = 0; i < 4; ++i) {
		lb_v.push_back(100);
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
	switch (m_obj){
		case MASS:
			f[0] = - x[11];
			break;
		case TIME:
			f[0] = (std::accumulate(x.begin() + 1,x.begin() + 8,0.) * ASTRO_DAY2YEAR);
			break;
		case MASS_TIME:
			f[0] = - x[11] / (std::accumulate(x.begin() + 1,x.begin() + 8,0.) * ASTRO_DAY2YEAR);
	}
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
		epoch start(std::accumulate(x.begin(),x.begin() + 2 * i + 1, 0.),epoch::MJD),
			end(std::accumulate(x.begin(),x.begin() + 2 * i + 2, 0.),epoch::MJD);
		// Set leg's start-end epochs.
		m_legs[i].set_t_i(start);
		m_legs[i].set_t_f(end);
		// Initial state.
		m_asteroids[i].eph(start,r,v);
		// First leg: Vinf is added to the state.
		if (i == 0) {
			for (int j = 0; j < 3; ++j) {
				v[j] += x[12 + j] * 1000;
			}
		}
		m_legs[i].set_x_i(sc_state(r,v,initial_mass));
		// Throttles.
		for (int j = 0; j < m_n_seg; ++j) {
			m_legs[i].set_throttles(j,get_nth_throttle(j,x.begin() + 15 + 3 * i * m_n_seg,start,end));
		}
		// Final state.
		m_asteroids[i + 1].eph(end,r,v);
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
		c[7*i+6] /=1500;
	}
	// Throttles constraints.
	for (int i = 0; i < 4; ++i) {
		m_legs[i].get_throttles_con(c.begin() + 28 + i * m_n_seg, c.begin() + 28 + (i + 1) * m_n_seg);
	}
	// Vinf constraint.
	c.back() = (x[12] * x[12] + x[13] * x[13] + x[14] * x[14] - 3.5 * 3.5) / (ASTRO_EARTH_VELOCITY * ASTRO_EARTH_VELOCITY) * 1000000;
}

/// Implementation of the sparsity structure: automated detection
//void gtoc_2::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
// {
//	//Numerical procedure
//	estimate_sparsity(lenG, iGfun, jGvar);
// }

std::string gtoc_2::get_name() const
{
	return "GTOC_2";
}

std::string gtoc_2::human_readable_extra() const {
	std::ostringstream oss;
	oss << "\n\tAsteroid Sequence:\t\t\t";
	for (int i =1 ; i< 5; ++i ) oss << m_asteroids[i].get_name() << " ";
	oss << std::endl;
	return oss.str();
}

/// A pretty description of the chromosome
/**
 * @return a string containing the human readable decision vector description
 */
std::string gtoc_2::pretty(const std::vector<double> &x) const
{
	std::ostringstream s;
	//We start by filling up the m_legs with the correct information
	constraint_vector c(this->get_c_dimension());
	this->compute_constraints_impl(c,x);
	s << "Final Mass: " << x[11] << " Kg" << std::endl;
	s << "Total Time: "  << std::accumulate(x.begin() + 1,x.begin() + 8,0.) * ASTRO_DAY2YEAR << " Years" << std::endl;
	s << "Objective Function: " << x[11] / ( std::accumulate(x.begin() + 1,x.begin() + 8,0.) * ASTRO_DAY2YEAR )  << " Kg/Years" << std::endl;
	for (int i=0; i<4;++i) s << this->m_legs[i] << std::endl;
	s << std::endl;

	return s.str();
}

}} // namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::gtoc_2)
