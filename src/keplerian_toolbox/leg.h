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

class state {
};

class dynamicalsystem {
	state s;
	double t;
	virtual void propagate(const double &DT) = 0;
	state forecast(const double &tnew) {
		dynamical_system new_s(*this);
		new_s.propagate(tnew - t);
		return new_s.s;
	}
}

class keplerianobject:dynamicalsystem {
	void propagate(const double &DT);
	state eph(const double &tnew) {
		return forecast(tnew);
	}
	
};

class endpoint {
	public:
		enum type {rndv, flyby}; 
		endpoint(type ty, const dynamicalsystem &ds):t(ty),s(ds.clone()) {}
		state eph(t) {
			return s.eph(t);
		}
		const dynamicalsystem &get_ds() const {
			return *s;
		}
	private:
		type t;
	        dynamicalsystem* 	s;
		bool	fixed;
}


class leg {
	leg(int n, const endpoint dep, const endpoint arr):departure(dep),arrival(arr),n_segments(n) {}
	endpoint		departure;
	endpoint		arrival;
	const int 	n_segments;
};

class earth_mars_lt_problem {
	earth_mars_lt_problem(int n, double thrust):l(n,endpoint(departure,keplerianobject(earth)),endpoint(rndv,keplerianobject(mars))),Tmax(thrust) {};
	
	double objfun(const std::vector<double> &x) {
		return main_objfun(x) + penalty(x);
	}

	private:
	double main_objfun(const std::vector<double> &x) const {
		double retval = 0;
		const int size = x.size();
		for (int i = 4; i < size - 4; i+=3) {
			retval += x[i];
		}
		return retval;
	}
	double penalty(const std::vector<double> &x) const {
		double retval = 0.;
		const double dt = x.back() / n_segments;
		for (int i = 0; i < n_segments; ++i) {
			if (norm(&x[4 + 3 * i]) / dt > Tmax / mass) {
				retval += std::abs(norm(&x[4 + 3 * i]) / dt - Tmax / mass);
			}
		}
		return retval + state_mismatch(x).norm() + max(0,(sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3])-Vmax));
	}
	kstate state_mismatch(const std::vector<double> &x) {
		const int n_seg_fwd = std::ceil((n_segments + 1) / 2.), n_seg_back = (n_segments + 1) / 2;
		const double dt = x.back() / n_segments;
		double cart[3];
		// We define two virtual spacecraft, one travelling forward from the first endpoint and a second one travelling
		// backward from the second endpoint
		keplerianobject spacecraft_fwd(l.departure.get_ds().eph(x[0]),x[0]);
		keplerianobject spacecraft_back(l.arrival.get_ds().eph(x[0] + x.back()),x[0] + x.back());
		// Forward propagation		
		// Launcher velocity increment
		spacecraft_fwd.kick(pol2cart(&x[1]));
		// First propagation
		spacecraft_fwd.propagate(dt/2.);
		// Other propagations
		for (int i = 0; i < n_seg_fwd; ++i) {
			spacecraft_fwd.kick(pol2cart(&x[4 + 3 * i]));
			if (i == n_seg_fwd - 1) {
				spacecraft_fwd.propagate(dt / 2.);
			} else {
				spacecraft_fwd.propagate(dt);
			}
		}
		// Backward propagation		
		// Arrival velocity increment
		spacecraft_back.kick(pol2cart(&x.back() - 3));
		// First propagation
		spacecraft_back.propagate(-dt/2.);
		// Other propagations
		for (int i = 0; i < n_seg_back; ++i) {
			spacecraft_back.kick(pol2cart(&x.back() - 6 - 3 * i));
			if (i == n_seg_back - 1) {
				spacecraft_back.propagate(-dt / 2.);
			} else {
				spacecraft_back.propagate(-dt);
			}
		}
		return spacecraft_back.s - spacecraft_fwd.s;
	}
	leg l;
	
};

double[3] pol2cart(const double[3] &r) {
	double retval[3];
	retval[0] = r[0] * std::cos(r[1]) * std::sin(r[2]);
	retval[1] = r[0] * std::cos(r[1]) * std::cos(r[2]);
	retval[2] = r[0] * std::sin(r[1]);
	return retval;
}

