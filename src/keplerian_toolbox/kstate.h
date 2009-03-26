#ifndef KEPLERIAN_TOOLBOX_KSTATE_H
#define KEPLERIAN_TOOLBOX_KSTATE_H

#include <cmath>
#include <iostream>

#include "s2t.h"
#include "t2s.h"
#include "types.h"

namespace keplerian_toolbox
{
	class kstate {
		public:
			kstate();
			template <class Vector>
			kstate(const Vector &, const Vector &, const double &, const double &);
			void propagate(const double &);
			// Getters and setters.
			const array_d3 &get_r() const;
			void set_r(const array_d3 &);
			template <class Vector>
			void set_r(const Vector &);
			const array_d3 &get_v() const;
			void set_v(const array_d3 &);
			template <class Vector>
			void set_v(const Vector &);
		private:
			friend std::ostream &operator<<(std::ostream &, const kstate &);
			array_d3	m_r;
			array_d3	m_v;
			double		m_t;
			double		m_mu;
	};

	inline kstate::kstate():m_t(0.),m_mu(1.)
	{
		m_r[0] = 1;
		m_r[1] = m_r[2] = 0;
		m_v[1] = 1;
		m_v[0] = m_v[2] = 0;
	}

	template <class Vector>
	inline kstate::kstate(const Vector &r, const Vector &v, const double &t, const double &mu):m_t(t),m_mu(mu)
	{
		if (mu <= 0) {
			// TODO: throw.
		}
		if (r.size() != 3 || v.size() != 3) {
			// TODO: throw.
		}
		if (r[0] == 0. && r[1] == 0. && r[2] == 0.) {
			// TODO: throw.
		}
		for (size_t i = 0; i < 3; ++i) {
			m_r[i] = r[i];
			m_v[i] = v[i];
		}
	}

	inline void kstate::propagate(const double &tf)
	{
		// Compute current absolute values for position and radial velocity.
		const double r0 = std::sqrt(m_r[0] * m_r[0] + m_r[1] * m_r[1] + m_r[2] * m_r[2]);
		const double vr0 = (m_v[0] * m_r[0] + m_v[1] * m_r[1] + m_v[2] * m_r[2]) / r0;
		const double v0 = std::sqrt(m_v[0] * m_v[0] + m_v[1] * m_v[1] + m_v[2] * m_v[2]);
		// Compute alpha.
		const double alpha = (2. * m_mu - v0 * v0 * r0) / r0;
		// Compute s by solving Kepler's equation.
		const double s = t2s(tf,r0,vr0,v0,m_t,m_mu);
		// F and G functions.
		double F, G, Fp, Gp;
		array_d3 old_r = m_r, old_v = m_v;
		// TODO: in the F&G funcs we must check that the new r does not go to zero.
		// Must decide a course of action: complete stop of the propagation (i.e., exception) or use some numerical trick to
		// set the radius to _almost_ zero?
		// Handle the parabolic case.
		if (alpha == 0) {
			const double s2 = s * s;
			// Calculate the F and G functions.
			F = 1. - m_mu / r0 * s2 / 2.;
			G = tf - m_t - m_mu * s * s2 / 6.;
			m_r[0] = F * old_r[0] + G * old_v[0];
			m_r[1] = F * old_r[1] + G * old_v[1];
			m_r[2] = F * old_r[2] + G * old_v[2];
			const double r = std::sqrt(m_r[0] * m_r[0] + m_r[1] * m_r[1] + m_r[2] * m_r[2]);
			// Now Fp and Gp.
			Fp = -m_mu * s / (r * r0);
			Gp = 1. - m_mu / r * s2 / 2.;
			m_v[0] = Fp * old_r[0] + Gp * old_v[0];
			m_v[1] = Fp * old_r[1] + Gp * old_v[1];
			m_v[2] = Fp * old_r[2] + Gp * old_v[2];
		} else {
			const double as2 = alpha * s * s;
			// Calculate the needed Stumpff functions.
			const double c0 = stumpff(0,as2), c1 = stumpff(1,as2);
			// Calculate the F and G functions.
			F = 1. - m_mu / r0 * (1. - c0) / alpha;
			G = tf - m_t - m_mu * s * (1. - c1) / alpha;
			m_r[0] = F * old_r[0] + G * old_v[0];
			m_r[1] = F * old_r[1] + G * old_v[1];
			m_r[2] = F * old_r[2] + G * old_v[2];
			const double r = std::sqrt(m_r[0] * m_r[0] + m_r[1] * m_r[1] + m_r[2] * m_r[2]);
			// Now Fp and Gp.
			Fp = -m_mu * s * c1 / (r * r0);
			Gp = 1. - m_mu / r * (1. - c0) / alpha;
			m_v[0] = Fp * old_r[0] + Gp * old_v[0];
			m_v[1] = Fp * old_r[1] + Gp * old_v[1];
			m_v[2] = Fp * old_r[2] + Gp * old_v[2];
		}
		// Last step: update the time.
		m_t = tf;
	}

	inline const array_d3 &kstate::get_r() const
	{
		return m_r;
	}

	inline void kstate::set_r(const array_d3 &r)
	{
		// TODO: check that r is not fucked up.
		m_r = r;
	}

	template <class Vector>
	inline void kstate::set_r(const Vector &r)
	{
		if (r.size() != 3) {
			// TODO: throw.
		}
		// TODO: check that r is not fucked up.
		m_r[0] = r[0];
		m_r[1] = r[1];
		m_r[2] = r[2];
	}

	inline const array_d3 &kstate::get_v() const
	{
		return m_v;
	}

	inline void kstate::set_v(const array_d3 &v)
	{
		m_v = v;
	}

	template <class Vector>
	inline void kstate::set_v(const Vector &v)
	{
		if (v.size() != 3) {
			// TODO: throw.
		}
		m_v[0] = v[0];
		m_v[1] = v[1];
		m_v[2] = v[2];
	}

	inline std::ostream &operator<<(std::ostream &os, const kstate &k)
	{
		os << "Position: [" << k.m_r[0] << ',' << k.m_r[1] << ',' << k.m_r[2] << "]\n";
		os << "Velocity: [" << k.m_v[0] << ',' << k.m_v[1] << ',' << k.m_v[2] << "]\n";
		os << "Time:     " << k.m_t << '\n';
		os << "Mu:       " << k.m_mu << '\n';
		return os;
	}
}

#endif
