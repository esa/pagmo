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

#ifndef PARTICLE_H
#define PARTICLE_H

#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <vector>

#include "coordinate_system.h"
#include "cartesian_cs.h"

// TODO: assert state size and positive mass, everywhere!
class particle {
		friend std::ostream &operator<<(std::ostream &, const particle &);
	public:
		/// Default constructor.
		/**
		 * Constructs a particle of unitary mass and zero velocity in the center of a cartesian coordinate system.
		 */
		particle():m_mass(1),m_state(6),m_cs(new cartesian_cs()) {}
		/// Constructor from mass, pointer to state vector and coordinate system.
		particle(const double &mass, const double *state, const coordinate_system &cs):
			m_mass(mass),m_state(state,state + 6),m_cs(cs.clone())
		{}
		/// Constructor from mass, state vector and coordinate system.
		particle(const double &mass, const std::vector<double> &state, const coordinate_system &cs):
			m_mass(mass),m_state(state),m_cs(cs.clone())
		{}
		/// Constructor from mass, state array and coordinate system.
		particle(const double &mass, const boost::array<double,6> &state, const coordinate_system &cs):
			m_mass(mass),m_state(state.begin(),state.end()),m_cs(cs.clone())
		{}
		/// Copy constructor.
		particle(const particle &p):m_mass(p.m_mass),m_state(p.m_state),m_cs(p.m_cs->clone()) {}
		/// Assignment operator.
		particle &operator=(const particle &p)
		{
			if (this != &p) {
				m_mass = p.m_mass;
				m_state = p.m_state;
				m_cs = p.m_cs->clone();
			}
			return *this;
		}
		/// Equality operator.
		bool operator==(const particle &p) const
		{
			return (m_mass == p.m_mass && m_state == p.m_state &&
				m_cs->operator==(*p.m_cs));
		}
		/// Set coordinate system.
		void set_cs(const coordinate_system &cs)
		{
			m_cs->to_cartesian(m_state);
			m_cs = cs.clone();
			m_cs->from_cartesian(m_state);
		}
		/// Get coordinate system.
		boost::shared_ptr<coordinate_system> get_cs() const
		{
			return m_cs->clone();
		}
		/// Get mass value.
		double get_mass() const
		{
			return m_mass;
		}
		/// Set mass value.
		void set_mass(const double &mass)
		{
			m_mass = mass;
		}
		/// Get const reference to state vector.
		const std::vector<double> &get_state() const
		{
			return m_state;
		}
		/// Set state vector value.
		void set_state(const std::vector<double> &state)
		{
			m_state = state;
		}
		/// Set state vector value.
		void set_state(const double *p_state)
		{
			for (size_t i = 0; i < 6; ++i) {
				m_state[i] = p_state[i];
			}
		}
		/// Set state vector value.
		void set_state(const boost::array<double,6> state)
		{
			for (size_t i = 0; i < 6; ++i) {
				m_state[i] = state[i];
			}
		}
	private:
		double					m_mass;
		std::vector<double>			m_state;
		boost::shared_ptr<coordinate_system>	m_cs;
};

std::ostream &operator<<(std::ostream &o, const particle &p)
{
	// TODO: assert state length of 6.
	o << std::scientific;
	o.precision(15);
	o << "mass: " << p.m_mass << '\n';
	o << "state: [";
	for (size_t i = 0; i < 5; ++i) {
		o << p.m_state[i] << ',';
	}
	o << p.m_state[5] << "]\n";
	o << "coordinate system: " << typeid(*(p.m_cs)).name() << '\n';
	return o;
}

#endif
