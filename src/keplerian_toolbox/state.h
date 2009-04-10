#ifndef KEPLERIAN_TOOLBOX_STATE_H
#define KEPLERIAN_TOOLBOX_STATE_H

#include <boost/array.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/static_assert.hpp>
#include <cmath>
#include <iostream>
#include <p_exceptions.h>

#include "types.h"

namespace keplerian_toolbox
{
	template <class T, int Size>
	class state {
			BOOST_STATIC_ASSERT(Size > 0);
		public:
			typedef T value_type;
			typedef size_t size_type;
			state() {
				for (size_type i = 0; i < Size; ++i) {
					m_array[i] = value_type(0);
				}
			}
			template <class Vector>
			state(const Vector &v) {
				if (v.size() != Size) {
					P_EX_THROW(value_error,"invalid vector size while constructing state");
				}
				for (size_type i = 0; i < Size; ++i) {
					m_array[i] = value_type(v[i]);
				}
			}
			size_type size() const {
				return Size;
			}
			const value_type &operator[](const size_type &n) const {
				return m_array[n];
			}
			value_type &operator[](const size_type &n) {
				return m_array[n];
			}
		private:
			boost::array<value_type,Size> m_array;
	};

	template <class T, int Size>
	inline std::ostream &operator<<(std::ostream &o, const state<T,Size> &s) {
		typedef typename state<T,Size>::size_type size_type;
		o << "State vector: [";
		for (size_type i = 0; i < Size; ++i) {
			o << s[i];
			if (i != Size - 1) {
				o << ' ';
			}
		}
		o << "]\n";
		return o;
	}

	template <class T>
	class cartesian_state;

	template <class T>
	class spherical_state;

	template <class T>
	class pv_state: public state<T,6> {
			typedef state<T,6> ancestor;
		public:
			pv_state():ancestor() {}
			template <class Vector>
			pv_state(const Vector &v):ancestor(v) {}
			template <class Vector>
			pv_state(const Vector &pos, const Vector &vel) {
				pv_size_check(pos);
				pv_size_check(vel);
				typedef typename ancestor::size_type size_type;
				for (size_type i = 0; i < 3; ++i) {
					(*this)[i] = pos[i];
					(*this)[i + 3] = vel[i];
				}
			}
			virtual ~pv_state() {}
			boost::array<T,3> get_position() const {
				boost::array<T,3> retval = {{(*this)[0], (*this)[1], (*this)[2]}};
				return retval;
			}
			boost::array<T,3> get_velocity() const {
				boost::array<T,3> retval = {{(*this)[3], (*this)[4], (*this)[5]}};
				return retval;
			}
			template <class Vector>
			void set_position(const Vector &p) {
				pv_size_check(p);
				(*this)[0] = p[0];
				(*this)[1] = p[1];
				(*this)[2] = p[2];
			}
			template <class Vector>
			void set_velocity(const Vector &v) {
				pv_size_check(v);
				(*this)[3] = v[0];
				(*this)[4] = v[1];
				(*this)[5] = v[2];
			}
			virtual cartesian_state<T> to_cartesian() const = 0;
			virtual spherical_state<T> to_spherical() const = 0;
		private:
			template <class Vector>
			static void pv_size_check(const Vector &v) {
				if (v.size() != 3) {
					P_EX_THROW(value_error,"invalid size for position/velocity vector");
				}
			}
	};

	template <class T>
	class cartesian_state: public pv_state<T> {
			typedef pv_state<T> ancestor;
		public:
			cartesian_state():ancestor() {}
			template <class Vector>
			cartesian_state(const Vector &v):ancestor(v) {}
			template <class Vector>
			cartesian_state(const Vector &pos, const Vector &vel):ancestor(pos,vel) {}
			virtual cartesian_state<T> to_cartesian() const {
				return *this;
			}
			virtual spherical_state<T> to_spherical() const;
	};

	template <class T>
	class spherical_state: public pv_state<T> {
			typedef pv_state<T> ancestor;
		public:
			spherical_state():ancestor() {}
			template <class Vector>
			spherical_state(const Vector &v):ancestor(v) {}
			template <class Vector>
			spherical_state(const Vector &pos, const Vector &vel):ancestor(pos,vel) {}
			virtual cartesian_state<T> to_cartesian() const;
			virtual spherical_state<T> to_spherical() const {
				return *this;
			}
	};

	template <class T>
	inline spherical_state<T> cartesian_state<T>::to_spherical() const {
		// Position.
		const T x = (*this)[0], y = (*this)[1], z = (*this)[2];
		const T r = std::sqrt(x * x + y * y + z * z);
		T phi, theta;
		// If r is zero, assign zero to the other guys by convention.
		if (r == 0) {
			phi = theta = T(0);
		} else {
			phi = std::atan2(y,x);
			theta = std::acos(z / r);
		}
		// Velocity.
		const T vx = (*this)[3], vy = (*this)[4], vz = (*this)[5];
		const T v = std::sqrt(vx * vx + vy * vy + vz * vz);
		T vphi, vtheta;
		if (v == 0) {
			vphi = vtheta = T(0);
		} else {
			vphi = std::atan2(vy,vx);
			vtheta = std::acos(vz / v);
		}
		boost::array<T,6> ret_array = {{r, phi, theta, v, vphi, vtheta}};
		return spherical_state<T>(ret_array);
	}

	template <class T>
	inline cartesian_state<T> spherical_state<T>::to_cartesian() const {
		// Position.
		const T r = (*this)[0], phi = (*this)[1], theta = (*this)[2], sin_theta = std::sin(theta),
			cos_theta = std::cos(theta), sin_phi = std::sin(phi), cos_phi = std::cos(phi);
		const T x = r * sin_theta * cos_phi, y = r * sin_theta * sin_phi, z = r * cos_theta;
		// Velocity.
		const T v = (*this)[0], vphi = (*this)[1], vtheta = (*this)[2], sin_vtheta = std::sin(vtheta),
			cos_vtheta = std::cos(vtheta), sin_vphi = std::sin(vphi), cos_vphi = std::cos(vphi);
		const T vx = v * sin_vtheta * cos_vphi, vy = v * sin_vtheta * sin_vphi, vz = v * cos_vtheta;
		boost::array<T,6> ret_array = {{x, y, z, vx, vy, vz}};
		return cartesian_state<T>(ret_array);
	}
}

#endif
