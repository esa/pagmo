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

#ifndef KEPLERIAN_TOOLBOX_STUMPFF_H
#define KEPLERIAN_TOOLBOX_STUMPFF_H

#include <boost/math/special_functions/factorials.hpp>
#include <cmath>

#define EPS (1e-15)

namespace keplerian_toolbox
{
	// Threshold for argument of Stumpff functions under which power series expansion is used.
	const double stumpff_series_thresh = .1;

	static inline int cs_phase(int n)
	{
		if (n & 1) {
			return -1;
		} else {
			return 1;
		}
	}

	static inline double stumpff_0(const double &x)
	{
		if (x >= 0) {
			return std::cos(std::sqrt(x));
		} else {
			return std::cosh(std::sqrt(-x));
		}
	}

	static inline double stumpff_1(const double &x)
	{
		if (x >= 0) {
			const double sqrt_x = std::sqrt(x);
			return std::sin(sqrt_x) / sqrt_x;
		} else {
			const double sqrt_x = std::sqrt(-x);
			return std::sinh(sqrt_x) / sqrt_x;
		}
	}

	static inline unsigned stumpff_series_limit(const double &x)
	{
		const double abs_x = std::abs(x);
		if (abs_x <= EPS) {
			return 1;
		}
		const double log_arg = std::log(1. - abs_x);
		if (log_arg == 0) {
			return 1;
		}
		return (unsigned)std::ceil((std::log(EPS) - std::log(-log_arg)) / std::log(abs_x));
	}

	static inline double stumpff_series(const unsigned &n, const double &x)
	{
		// This is the limit of the power series expansion.
		const unsigned int p = stumpff_series_limit(x);
		double retval = 0;
		double x_pow_i = 1.;
		double n_2i_fact = boost::math::factorial<double>(n);
		// Here the boundary is <= because in the original definition the summation starts from
		// 0 to p included.
		for (unsigned i = 0; i <= p; ++i) {
			retval += cs_phase(i) * (x_pow_i / n_2i_fact);
			// Next values for x_pow_i and the dividing factorial.
			x_pow_i *= x;
			const double tmp = n + (i + 1.) * 2.;
			n_2i_fact *= tmp * (tmp - 1.);
		}
		return retval;
	}

	inline double stumpff(int n, const double &x)
	{
		const unsigned n_ = n;
		const double abs_x = std::abs(x);
		if (abs_x == 0) {
			return (1. / boost::math::factorial<double>(n));
		}
		if (abs_x < stumpff_series_thresh) {
			return stumpff_series(n_,x);
		}
		if (n_ == 0) {
			return stumpff_0(x);
		}
		if (n_ == 1) {
			return stumpff_1(x);
		}
		const unsigned start_n = n_ & 1;
		unsigned cur_n = start_n;
		// The initial factorial will be either 0! == 1 or 1! == 1.
		double cur_fact = 1;
		double cur_stumpff;
		if (cur_n == 0) {
			cur_stumpff = stumpff_0(x);
		} else {
			cur_stumpff = stumpff_1(x);
		}
		for (unsigned i = cur_n; i != n_; i += 2) {
			cur_stumpff = (1. / cur_fact - cur_stumpff) / x;
			cur_fact *= (i + 2.) * (i + 1.);
		}
		return cur_stumpff;
	}

	inline double stumpff_0_p(const double &x)
	{
		const double abs_x = std::abs(x);
		double retval = -.5;
		if (abs_x == 0) {
			return retval;
		}
		if (abs_x < stumpff_series_thresh) {
			const unsigned int p = stumpff_series_limit(x);
			double x_pow_i = x;
			double fact = 24;
			for (unsigned i = 1; i <= p; ++i) {
				retval += cs_phase(i + 1) * (i + 1.) * (x_pow_i / fact);
				x_pow_i *= x;
				const double tmp = (i + 2.) * 2.;
				fact *= tmp * (tmp - 1.);
			}
			return retval;
		}
		if (x >= 0) {
			const double sqrt_x = std::sqrt(x);
			return -std::sin(sqrt_x) / (2. * sqrt_x);
		} else {
			const double sqrt_x = std::sqrt(-x);
			return -std::sinh(sqrt_x) / (2. * sqrt_x);
		}
	}
}

#undef EPS

#endif
