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

#ifndef KEPLERIAN_TOOLBOX_OE2STATE_H
#define KEPLERIAN_TOOLBOX_OE2STATE_H

#include <cmath>

#include "types.h"

namespace keplerian_toolbox
{
	inline void oe2state(const array_d6 &oe, const double &t, const double &mu, array_d3 &r, array_d3 &v)
	{
		// Cache values from vector of orbital elements.
		const double a = oe[0], e = oe[1], I = oe[2], o = oe[3], O = oe[4], M0 = oe[5];
		// TODO: check for bogus values of a (0), eccentricity (<0), I (-pi/2,pi/2, I suppose) and reduce o and O to
		// be in [0,2pi[. Also, check for mu and r as usual. Also, check that a and e are consistent...
		// Maybe a general checker function for orbital elements can be used here and in kstate. This function could
		// take the input oe vector and modify it if needed (to put the angles in the right ranges, etc.) and then feed it here.
		if (e >= 1) {
			std::cout << "Only elliptical orbits are supported at the moment.\n";
			std::exit(1);
		}
		const double n = std::sqrt(mu / (a * a * a));
	}
}

#endif
