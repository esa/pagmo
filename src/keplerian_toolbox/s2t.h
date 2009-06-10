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

#ifndef KEPLERIAN_TOOLBOX_S2T_H
#define KEPLERIAN_TOOLBOX_S2T_H

#include <cmath>

#include "../p_exceptions.h"
#include "stumpff.h"

namespace keplerian_toolbox
{
	inline double s2t(const double &s, const double &r0, const double &vr0, const double &v0, const double &t0, const double &mu)
	{
		if (mu <= 0) {
			P_EX_THROW(value_error,"mu must be strictly positive");
		}
		if (r0 <= 0) {
			P_EX_THROW(value_error,"initial radius must be strictly positive");
		}
		// Cache values.
		const double s2 = s * s, s3 = s2 * s, as2 = (s2 * (2. * mu - r0 * v0 * v0)) / r0;
		// Evaluate first 2 Stumpff functions.
		const double c0 = stumpff(0,as2), c1 = stumpff(1,as2);
		// Evaluate c2 and c3. Calculate from scratch if as2 is small, use recursive relations otherwise.
		double c2_, c3_;
		if (std::abs(as2) < stumpff_series_thresh) {
			c2_ = stumpff(2,as2);
			c3_ = stumpff(3,as2);
		} else {
			c2_ = (1. - c0) / as2;
			c3_ = (1. - c1) / as2;
		}
		// Recast to const for safety.
		const double c2 = c2_, c3 = c3_;
		return (t0 + r0 * s * c1 + r0 * vr0 * s2 * c2 + mu * s3 * c3);
	}
}

#endif
