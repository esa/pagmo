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

// 05/2009: initial version by Dario Izzo and Francesco Biscani.

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include "../../AstroToolbox/Pl_Eph_An.h"
#include "../../AstroToolbox/propagateKEP.h"
#include "../../constants.h"
#include "base.h"
#include "earth_mars_lt.h"

#define R 149597870.66
#define V std::sqrt(1.32712428e+11 / 149597870.66)
#define T (R / V)
#define A (R / (T * T))
#define F (1 * A)

namespace pagmo
{
namespace problem {

// NOTE: parameters must be dimensional:
// - thrust in Newton
// - mass in kilograms
earth_mars_lt::earth_mars_lt(int n_, const double &M_, const double &thrust_, const double &Isp_):
		base(n_ * 3 + 5),n(n_),M(M_),thrust((thrust_ / 1000) / F),Isp(Isp_)
{
// std::cout << "thrust:\t" << thrust << '\n';
// std::cout << R << '\n';
// std::cout << V << '\n';
// std::cout << T << '\n';
// std::cout << A << '\n';
// std::cout << F << '\n';
// std::cout << M << '\n';
// std::cout << "-------" << '\n';
	// NOTE: sphere picking coord.
	// Dep. date.
	LB[0] = 0;
	UB[0] = 5000;
	// Vinf at departure.
	LB[1] = 0;
	UB[1] = 3;
	LB[2] = 0;
	UB[2] = 1;
	LB[3] = 0;
	UB[3] = 1;
	for (int i = 0; i < n; ++i) {
		// Segments' deltaVs.
		LB[3 * i + 4] = 0;
		UB[3 * i + 4] = 1;
		LB[3 * i + 5] = 0;
		UB[3 * i + 5] = 1;
		LB[3 * i + 6] = 0;
		UB[3 * i + 6] = 1;
	}
	// Transfer time.
	LB.back() = 0;
	UB.back() = 1000;
}

void earth_mars_lt::human_readable(const std::vector<double> &x) const
{
	using namespace std;
	cout << left << "*********************** " << endl;
	cout << "Trajectory description: " << endl;
	cout << "*********************** " << endl << endl;
	cout << setw(40) << "Number of segments:" << n << endl <<
	     setw(40) << "Max thrust:" << thrust * F * 1000 << " N" << endl <<
	     setw(40) << "Max DV per segment:" << (thrust * F / M * x.back() * 86400. / n) << " Km/s" << endl <<
	     setw(40) << "Initial mass:" << M << " Kg" << endl <<
	     setw(40) << "Departure date:" << x[0] << " MJD2000" << endl <<
	     setw(40) << "Vinf at departure:" << x[1] << " Km/s" << endl << endl;

	cout << setw(25) << "DV" << setw(15) << right << "vx" << setw(15) << "vy" << setw(15) << "vz" << left << endl;
	double tmp_velocity[3];
	double dt = (x.back() / n) * 86400 / T;
	for (size_t i = 0; i < (size_t)n; ++i) {
		ruv2cart(tmp_velocity,&x[3 * i + 4]);
		cout << setw(25) << x[3 * i + 4] * thrust / M * dt * V << right;
		cout << setw(15) << tmp_velocity[0] * thrust / M * dt * V << setw(15) << tmp_velocity[1] * thrust / M * dt * V << setw(15) << tmp_velocity[2] * thrust / M * dt * V << left << endl;
	}


	double r_fwd[3], v_fwd[3], r_back[3], v_back[3],Dr,Dv;
	state_mismatch(x,r_fwd,v_fwd,r_back,v_back);
	Dr = (r_fwd[0] - r_back[0])*(r_fwd[0] - r_back[0]) +(r_fwd[1] - r_back[1])*(r_fwd[1] - r_back[1]) + (r_fwd[2] - r_back[2])*(r_fwd[2] - r_back[2]);
	Dv = (v_fwd[0] - v_back[0])*(v_fwd[0] - v_back[0]) +(v_fwd[1] - v_back[1])*(v_fwd[1] - v_back[1]) + (v_fwd[2] - v_back[2])*(v_fwd[2] - v_back[2]);
	Dr = sqrt(Dr);
	Dv = sqrt(Dv);
	cout << endl << setw(25) << "Pos mismatch (km):" << right << setw(15) << Dr*R << left;// << setw(15) << (r_fwd[0] - r_back[0]) << setw(15) << (r_fwd[1] - r_back[1]) << setw(15) << (r_fwd[2] - r_back[2]) << left << endl;
	cout << endl << setw(25) << "Vel mismatch (km/sec):" << right << setw(15) << Dv*V << left << endl << endl;// << setw(15) << (v_fwd[0] - v_back[0]) << setw(15) << (v_fwd[1] - v_back[1]) << setw(15) << (v_fwd[2] - v_back[2]) << left << endl << endl;

	cout << setw(40) << "Time of flight:" << x.back() << " days" << endl;
	cout << setw(40) << "Total DV:" << main_objfun(x) * V << " Km/s" << endl;
	cout << setw(40) << "Final mass:" << M * std::exp(-main_objfun(x) * V * 1000. / (9.80665 * Isp)) << " Kg" << endl;
}

double earth_mars_lt::objfun_(const std::vector<double> &x) const
{
	double r_fwd[3], v_fwd[3], r_back[3], v_back[3];
	state_mismatch(x,r_fwd,v_fwd,r_back,v_back);
	const double s_mismatch = std::sqrt((r_back[0] - r_fwd[0]) * (r_back[0] - r_fwd[0]) + (r_back[1] - r_fwd[1]) * (r_back[1] - r_fwd[1]) +
	                                    (v_back[0] - v_fwd[0]) * (v_back[0] - v_fwd[0]) + (v_back[1] - v_fwd[1]) * (v_back[1] - v_fwd[1]) +
	                                    (v_back[2] - v_fwd[2]) * (v_back[2] - v_fwd[2]));
//std::cout << s_mismatch << '\n';
	return main_objfun(x) + 1000 * s_mismatch;
}

double earth_mars_lt::main_objfun(const std::vector<double> &x) const
{
	double retval = 0;
	double dt = (x.back() / n) * 86400 / T;
	for (int i = 0; i < n; ++i) {
		retval += x[3 * i + 4] * thrust / M * dt ;
	}
	return retval;
}

void earth_mars_lt::state_mismatch(const std::vector<double> &x, double *r_fwd, double *v_fwd, double *r_back, double *v_back) const
{
// std::cout << "max DV:\t" << thrust / M * (x.back() / n) * 86400 / T << '\n';
	const int n_seg_fwd = (n + 1) / 2, n_seg_back = n / 2;
// std::cout << "n segments:\t" << n_seg_fwd << ',' << n_seg_back << '\n';
	const double dt = (x.back() / n) * 86400 / T;
// std::cout << "dt_initial:\t" << dt << '\n';
	// We define two virtual spacecraft, one travelling forward from the first endpoint and a second one travelling
	// backward from the second endpoint.
	//double r_fwd[3], v_fwd[3], r_back[3], v_back[3];
	earth_eph(x[0],r_fwd,v_fwd);
// std::cout << "rs_fwd:\t" << r_fwd[0] << ',' << r_fwd[1] << ',' << r_fwd[2] << '\n';
// std::cout << "vs_fwd:\t" << v_fwd[0] << ',' << v_fwd[1] << ',' << v_fwd[2] << '\n';
	mars_eph(x[0] + x.back(),r_back,v_back);
	// Forward propagation.
	// Launcher velocity increment
	kick(v_fwd,&x[1]);
// std::cout << "rs_fwd:\t" << r_fwd[0] << ',' << r_fwd[1] << ',' << r_fwd[2] << '\n';
// std::cout << "vs_fwd:\t" << v_fwd[0] << ',' << v_fwd[1] << ',' << v_fwd[2] << '\n';
	// First propagation
	propagate(r_fwd,v_fwd,dt / 2.);
// std::cout << "rs_fwd:\t" << r_fwd[0] << ',' << r_fwd[1] << ',' << r_fwd[2] << '\n';
// std::cout << "vs_fwd:\t" << v_fwd[0] << ',' << v_fwd[1] << ',' << v_fwd[2] << '\n';
	// Other propagations
	for (int i = 0; i < n_seg_fwd; ++i) {
		punch(v_fwd,&x[3 * i + 4],thrust / M,dt);
// std::cout << "rs_fwd:\t" << r_fwd[0] << ',' << r_fwd[1] << ',' << r_fwd[2] << '\n';
// std::cout << "vs_fwd:\t" << v_fwd[0] << ',' << v_fwd[1] << ',' << v_fwd[2] << '\n';
		if (i == n_seg_fwd - 1) {
			propagate(r_fwd,v_fwd,dt / 2.);
// std::cout << "rs_fwd:\t" << r_fwd[0] << ',' << r_fwd[1] << ',' << r_fwd[2] << '\n';
// std::cout << "vs_fwd:\t" << v_fwd[0] << ',' << v_fwd[1] << ',' << v_fwd[2] << '\n';
		} else {
			propagate(r_fwd,v_fwd,dt);
// std::cout << "rs_fwd:\t" << r_fwd[0] << ',' << r_fwd[1] << ',' << r_fwd[2] << '\n';
// std::cout << "vs_fwd:\t" << v_fwd[0] << ',' << v_fwd[1] << ',' << v_fwd[2] << '\n';
		}
	}
	// Backward propagation.
	// First propagation
	back_propagate(r_back,v_back,dt/2.);
// std::cout << "rs_back:\t" << r_back[0] << ',' << r_back[1] << ',' << r_back[2] << '\n';
// std::cout << "vs_back:\t" << v_back[0] << ',' << v_back[1] << ',' << v_back[2] << '\n';
	// Other propagations
	for (int i = 0; i < n_seg_back; ++i) {
		back_punch(v_back,&x.back() - 3 - 3 * i,thrust / M,dt);
// std::cout << "rs_back:\t" << r_back[0] << ',' << r_back[1] << ',' << r_back[2] << '\n';
// std::cout << "vs_back:\t" << v_back[0] << ',' << v_back[1] << ',' << v_back[2] << '\n';
		if (i == n_seg_back - 1) {
			back_propagate(r_back,v_back,dt / 2.);
// std::cout << "rs_back:\t" << r_back[0] << ',' << r_back[1] << ',' << r_back[2] << '\n';
// std::cout << "vs_back:\t" << v_back[0] << ',' << v_back[1] << ',' << v_back[2] << '\n';
		} else {
			back_propagate(r_back,v_back,dt);
// std::cout << "rs_back:\t" << r_back[0] << ',' << r_back[1] << ',' << r_back[2] << '\n';
// std::cout << "vs_back:\t" << v_back[0] << ',' << v_back[1] << ',' << v_back[2] << '\n';
		}
	}
// 			return std::sqrt((r_back[0] - r_fwd[0]) * (r_back[0] - r_fwd[0]) + (r_back[1] - r_fwd[1]) * (r_back[1] - r_fwd[1]) +
// 				(r_back[2] - r_fwd[2]) * (r_back[2] - r_fwd[2])) +
// 				std::sqrt((v_back[0] - v_fwd[0]) * (v_back[0] - v_fwd[0]) + (v_back[1] - v_fwd[1]) * (v_back[1] - v_fwd[1]) +
// 					(v_back[2] - v_fwd[2]) * (v_back[2] - v_fwd[2]));
}

void earth_mars_lt::ruv2cart(double *output, const double *input)
{
	double theta, phi, r = input[0];
	theta = 2 * M_PI * input[1];
	phi = std::acos(2 * input[2] - 1);
	output[0] = r * std::cos(theta) * std::sin(phi);
	output[1] = r * std::sin(theta) * std::sin(phi);
	output[2] = r * std::cos(phi);
}

void earth_mars_lt::earth_eph(const double &mjd2000, double *position, double *velocity)
{
	// Non-dimensional axis!
	const double keplerian[6] = {149597891.091458 / R , 0.0167090220204814, 0, 0, -257.059521691888, -2.47457303222881};
	// Epoch must be in mjd.
	const double epoch = 51544;
	Custom_Eph(mjd2000 + 2.451544500000000e+06, epoch, keplerian, position, velocity);
	position[0] /= R;
	position[1] /= R;
	position[2] /= R;
	velocity[0] /= V;
	velocity[1] /= V;
	velocity[2] /= V;
// std::cout << "Earth pos:\t" << position[0] << ',' << position[1] << ',' << position[2] << '\n';
// std::cout << "Earth vel:\t" << velocity[0] << ',' << velocity[1] << ',' << velocity[2] << '\n';
}

void earth_mars_lt::mars_eph(const double &mjd2000, double *position, double *velocity)
{
	// Non-dimensional axis!
	const double keplerian[6] = {227940515.511383 / R, 0.0934047954172235, 1.84967094444443, 49.5574266111111, -73.4983588723774,19.3881291190116};
	// Epoch must be in mjd.
	const double epoch = 51544;
	Custom_Eph(mjd2000 + 2.451544500000000e+06, epoch, keplerian, position,velocity);
	position[0] /= R;
	position[1] /= R;
	position[2] /= R;
	velocity[0] /= V;
	velocity[1] /= V;
	velocity[2] /= V;
// std::cout << "Mars pos:\t" << position[0] << ',' << position[1] << ',' << position[2] << '\n';
// std::cout << "Mars vel:\t" << velocity[0] << ',' << velocity[1] << ',' << velocity[2] << '\n';
}

void earth_mars_lt::kick(double *output, const double *input)
{
	double input_cart[3], input_copy[3];
	input_copy[0] = input[0] / V;
	input_copy[1] = input[1];
	input_copy[2] = input[2];
	ruv2cart(input_cart,input_copy);
	output[0] += input_cart[0];
	output[1] += input_cart[1];
	output[2] += input_cart[2];
}

void earth_mars_lt::punch(double *output, const double *input, const double &Tmax_M, const double &dt)
{
	double input_cart[3], input_copy[3];
	input_copy[0] = input[0] * dt * Tmax_M;
	input_copy[1] = input[1];
	input_copy[2] = input[2];
	ruv2cart(input_cart,input_copy);
	output[0] += input_cart[0];
	output[1] += input_cart[1];
	output[2] += input_cart[2];
}

void earth_mars_lt::back_punch(double *output, const double *input, const double &Tmax_M, const double &dt)
{
	double input_cart[3], input_copy[3];
	input_copy[0] = input[0] * dt * Tmax_M;
	input_copy[1] = input[1];
	input_copy[2] = input[2];
	ruv2cart(input_cart,input_copy);
	output[0] -= input_cart[0];
	output[1] -= input_cart[1];
	output[2] -= input_cart[2];
}

void earth_mars_lt::propagate(double *r, double *v, const double &t)
{
// std::cout << "dt:\t" << t << '\n';
	double dummy_r[3], dummy_v[3];
	propagateKEP(r,v,t,1,dummy_r,dummy_v);
	r[0] = dummy_r[0];
	r[1] = dummy_r[1];
	r[2] = dummy_r[2];
	v[0] = dummy_v[0];
	v[1] = dummy_v[1];
	v[2] = dummy_v[2];
}

void earth_mars_lt::back_propagate(double *r, double *v, const double &t)
{
// std::cout << "dt:\t" << t << '\n';
	double dummy_r[3], dummy_v[3], back_v[3];
	back_v[0] = -v[0];
	back_v[1] = -v[1];
	back_v[2] = -v[2];
	propagateKEP(r,back_v,t,1,dummy_r,dummy_v);
	r[0] = dummy_r[0];
	r[1] = dummy_r[1];
	r[2] = dummy_r[2];
	v[0] = -dummy_v[0];
	v[1] = -dummy_v[1];
	v[2] = -dummy_v[2];
}

}
}
