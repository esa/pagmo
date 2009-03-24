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
