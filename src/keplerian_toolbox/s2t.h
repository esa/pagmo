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
