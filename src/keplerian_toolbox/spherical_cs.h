#ifndef SPHERICAL_CS_H
#define SPHERICAL_CS_H

#include <boost/shared_ptr.hpp>
#include <cmath>
#include <vector>

#include "coordinate_system.h"

struct spherical_cs: public coordinate_system {
	void to_cartesian(std::vector<double> &s) const
	{
		// TODO: check p and v size.
		// Position.
		const double r = s[0], phi = s[1], theta = s[2], sin_theta = std::sin(theta),
			cos_theta = std::cos(theta), sin_phi = std::sin(phi), cos_phi = std::cos(phi);
		const double x = r * sin_theta * cos_phi, y = r * sin_theta * sin_phi, z = r * cos_theta;
		// Velocity.
		const double v_abs = s[3], vphi = s[4], vtheta = s[5], sin_vtheta = std::sin(vtheta),
			cos_vtheta = std::cos(vtheta), sin_vphi = std::sin(vphi), cos_vphi = std::cos(vphi);
		const double vx = v_abs * sin_vtheta * cos_vphi, vy = v_abs * sin_vtheta * sin_vphi, vz = v_abs * cos_vtheta;
		s[0] = x;
		s[1] = y;
		s[2] = z;
		s[3] = vx;
		s[4] = vy;
		s[5] = vz;
	}
	void from_cartesian(std::vector<double> &s) const
	{
		// TODO: check p and v size.
		// Position.
		const double x = s[0], y = s[1], z = s[2];
		const double r = std::sqrt(x * x + y * y + z * z);
		double phi, theta;
		// If r is zero, assign zero to the other guys by convention.
		if (r == 0) {
			phi = theta = 0;
		} else {
			phi = std::atan2(y,x);
			theta = std::acos(z / r);
		}
		// Velocity.
		const double vx = s[3], vy = s[4], vz = s[5];
		const double v_abs = std::sqrt(vx * vx + vy * vy + vz * vz);
		double vphi, vtheta;
		if (v_abs == 0) {
			vphi = vtheta = 0;
		} else {
			vphi = std::atan2(vy,vx);
			vtheta = std::acos(vz / v_abs);
		}
		s[0] = r;
		s[1] = phi;
		s[2] = theta;
		s[3] = v_abs;
		s[4] = vphi;
		s[5] = vtheta;
	}
	boost::shared_ptr<coordinate_system> clone() const
	{
		return boost::shared_ptr<coordinate_system>(new spherical_cs());
	}
};

#endif
