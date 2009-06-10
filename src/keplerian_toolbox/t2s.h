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

#ifndef KEPLERIAN_TOOLBOX_T2S_H
#define KEPLERIAN_TOOLBOX_T2S_H

#include <cmath>
#include <iostream>

#include "../p_exceptions.h"
#include "stumpff.h"

namespace keplerian_toolbox
{
	// Solve Kepler's equation in universal variable formulation using Laguerre's method.
	inline double t2s(const double &t, const double &r0, const double &vr0, const double &v0, const double &t0, const double &mu)
	{
		if (mu <= 0) {
			P_EX_THROW(value_error,"mu must be strictly positive");
		}
		if (r0 <= 0) {
			P_EX_THROW(value_error,"initial radius must be strictly positive");
		}
		// Configuration parameters.
		static const double n = 5., tol = 1E-13;
		static const unsigned max_it_n = 100;
		// Cache values.
		const double r0vr0 = r0 * vr0, alpha = (2. * mu - v0 * v0 * r0) / r0, delta_t = t0 - t, r0vr0_alpha = r0vr0 / alpha, a = mu / alpha, r0_m_a = r0 - a;
		// Initial guess for return value.
		double s = 0.;
		// In case of elliptical orbits, we can make an educated guess.
		if (alpha > 0) {
			s = (alpha * (t - t0) - r0vr0) / mu;
		}
		unsigned int i = 0;
		while (true) {
			const double s2 = s * s, as2 = alpha * s2;
			// Calculate f(s), fp(s), fpp(s).
			double f, fp, fpp;
			if (alpha == 0) {
				f = delta_t + r0 * s + r0vr0 * s2 / 2. + s2 * s * mu / 6.;
				fp = r0 + r0vr0 * s + mu * s2 / 2.;
				fpp = r0vr0 + mu * s;
			} else {
				// Calculate c0, c1, c0_p.
				const double c0 = stumpff(0,as2), c1 = stumpff(1,as2), c0_p = stumpff_0_p(as2);
				f = delta_t + r0vr0_alpha + a * s - c0 * r0vr0_alpha + s * c1 * r0_m_a;
				fp = a + c0 * r0_m_a + c1 * s * r0vr0;
				fpp = r0 * c0_p + r0vr0 * c0 + mu * s * c1;
			}
			// If f is small enough, stop.
			if (std::abs(f) <= tol) {
				break;
			}
			double cur_n = n;
			const double G = fp / f, G2 = G * G, H = G2 - fpp / f;
			double root = std::sqrt(std::abs((cur_n - 1.) * (cur_n * H - G2)));
			// Here we change root using another polynomial order n in order to avoid a division
			// by zero when G is also zero.
			while (G == 0 && root == 0) {
				cur_n += 1.;
				root = std::sqrt(std::abs((cur_n - 1.) * (cur_n * H - G2)));
			}
			double div = G - root, alt_div = G + root;
			if (std::abs(alt_div) > std::abs(div)) {
				div = alt_div;
			}
			const double diff = cur_n / div;
			s -= diff;
			++i;
			if (std::abs(diff) < tol) {
				break;
			}
			if (i > max_it_n) {
				std::cout << "Maximum number of iterations in t2s exceeded, exiting.\n";
				std::cout << "Current diff is: " << diff << '\n';
				std::cout << "Current f is: " << f << '\n';
				break;
			}
		}
		return s;
	}
}

#endif
