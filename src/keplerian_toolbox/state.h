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

#ifndef KEPLERIAN_TOOLBOX_STATE_H
#define KEPLERIAN_TOOLBOX_STATE_H

#include <boost/array.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/static_assert.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <cmath>
#include <iostream>
#include <vector>

#include "../p_exceptions.h"
#include "types.h"

namespace keplerian_toolbox
{
	template <class T, int Size>
	class state {
			BOOST_STATIC_ASSERT(Size > 0);
		public:
			typedef T value_type;
			typedef size_t size_type;
			state()
			{
				for (size_type i = 0; i < Size; ++i) {
					m_array[i] = value_type(0);
				}
			}
			template <class Vector>
			state(const Vector &v)
			{
				if (v.size() != Size) {
					P_EX_THROW(value_error,"invalid vector size while constructing state");
				}
				for (size_type i = 0; i < Size; ++i) {
					m_array[i] = value_type(v[i]);
				}
			}
			size_type size() const
			{
				return Size;
			}
			const value_type &operator[](const size_type &n) const
			{
				return m_array[n];
			}
		protected:
			boost::array<value_type,Size> m_array;
	};

	template <class T, int Size>
	inline std::ostream &operator<<(std::ostream &o, const state<T,Size> &s)
	{
		typedef typename state<T,Size>::size_type size_type;
		o << std::scientific;
		o.precision(15);
		o << "State vector: [";
		for (size_type i = 0; i < Size; ++i)
		{
			o << s[i];
			if (i != Size - 1) {
				o << ' ';
			}
		}
		o << "]\n";
		return o;
	}

	template <class T>
	struct base_coordinate_system {
		virtual void to_cartesian(boost::array<T,6> &) const {}
		virtual void from_cartesian(boost::array<T,6> &) const {}
		virtual boost::shared_ptr<base_coordinate_system> clone() const = 0;
		virtual ~base_coordinate_system() {}
	};

	template <class T>
	struct cartesian_coordinate_system: public base_coordinate_system<T> {
		boost::shared_ptr<base_coordinate_system<T> > clone() const
		{
			return boost::shared_ptr<base_coordinate_system<T> >(new cartesian_coordinate_system());
		}
	};

	template <class T>
	struct spherical_coordinate_system: public base_coordinate_system<T> {

	};

	template <class T>
	class pv_state: public state<T,6> {
			typedef state<T,6> ancestor;
		public:
			pv_state():ancestor(),m_cs(new cartesian_coordinate_system<T>()) {}
			template <class Vector>
			pv_state(const Vector &v):ancestor(v),m_cs(new cartesian_coordinate_system<T>()) {}
			template <class Vector>
			pv_state(const Vector &pos, const Vector &vel):ancestor(),m_cs(new cartesian_coordinate_system<T>())
			{
				pv_size_check(pos);
				pv_size_check(vel);
				typedef typename ancestor::size_type size_type;
				for (size_type i = 0; i < 3; ++i) {
					this->m_array[i] = pos[i];
					this->m_array[i + 3] = vel[i];
				}
			}
			pv_state(const pv_state &pvs):ancestor(pvs),m_cs(pvs.m_cs->clone()) {}
			pv_state &operator=(const pv_state &pvs)
			{
				if (this != &pvs) {
					ancestor::operator=(pvs);
					m_cs = pvs.m_cs->clone();
				}
				return *this;
			}
			boost::array<T,3> get_position() const
			{
				boost::array<T,3> retval = {{(*this)[0], (*this)[1], (*this)[2]}};
				return retval;
			}
			boost::array<T,3> get_velocity() const
			{
				boost::array<T,3> retval = {{(*this)[3], (*this)[4], (*this)[5]}};
				return retval;
			}
			boost::shared_ptr<base_coordinate_system<T> > get_coordinate_system() const
			{
				return m_cs->clone();
			}
			pv_state &set_coordinate_system(const base_coordinate_system<T> &cs)
			{
				m_cs->to_cartesian(this->m_array);
				m_cs = cs.clone();
				m_cs->from_cartesian(this->m_array);
				return *this;
			}
			template <class Vector>
			void set_position(const Vector &p)
			{
				pv_size_check(p);
				this->m_array[0] = p[0];
				this->m_array[1] = p[1];
				this->m_array[2] = p[2];
			}
			template <class Vector>
			void set_velocity(const Vector &v)
			{
				pv_size_check(v);
				this->m_array[3] = v[0];
				this->m_array[4] = v[1];
				this->m_array[5] = v[2];
			}
		private:
			template <class Vector>
			static void pv_size_check(const Vector &v)
			{
				if (v.size() != 3) {
					P_EX_THROW(value_error,"invalid size for position/velocity vector");
				}
			}
			boost::shared_ptr<base_coordinate_system<T> > m_cs;
	};

	template <class T>
	struct base_orbit_propagator {
		virtual void propagate(std::vector<pv_state<T> > &) const {}
		virtual bool verify(const std::vector<pv_state<T> > &) const
		{
			return true;
		}
		virtual boost::shared_ptr<base_orbit_propagator> clone() const = 0;
		virtual ~base_orbit_propagator() {}
	};

	template <class T>
	struct null_orbit_propagator: public base_orbit_propagator<T> {
		boost::shared_ptr<base_orbit_propagator<T> > clone() const
		{
			return boost::shared_ptr<base_orbit_propagator<T> >(new null_orbit_propagator());
		}
	};

	template <class T>
	class dynamical_system {
		public:
			dynamical_system():m_op(new null_orbit_propagator<T>()),m_states(),m_time(0)
			{
				if (!m_op->verify(m_states)) {
					
				}
			}
			dynamical_system(const dynamical_system &b):m_op(b.m_op->clone()),m_states(b.m_states),m_time(b.m_time) {}
			dynamical_system &operator=(const dynamical_system &b)
			{
				if (this != &b) {
					m_op = b.m_op->clone();
					m_states = b.m_states;
					m_time = b.m_time;
				}
				return *this;
			}
			size_t size() const
			{
				return m_states.size();
			}
		private:
			boost::shared_ptr<base_orbit_propagator<T> >	m_op;
			std::vector<pv_state<T>	>			m_states;
			double						m_time;
	};
}

#endif
