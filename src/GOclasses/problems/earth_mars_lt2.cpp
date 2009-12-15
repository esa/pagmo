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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include "../../AstroToolbox/Pl_Eph_An.h"
#include "../../AstroToolbox/taylor_fixedthrust.h"
#include "../../constants.h"
#include "../../exceptions.h"
#include "base.h"
#include "earth_mars_lt2.h"

#define R 149597870.66
#define V std::sqrt(1.32712428e+11 / 149597870.66)
#define T (R / V)
#define A (R / (T * T))
#define F (1 * A)

namespace pagmo
{
namespace problem {

int* param;

earth_mars_lt2::earth_mars_lt2(int n_, const double &M_, const double &thrust_, const double &Isp_):
		base(n_ * 3 + 5),n(n_),M(M_),thrust((thrust_ / 1000) / F),Isp(Isp_),eps(1e-5)
{

	// Launch date in MJD200
	LB[0] = 0;
	UB[0] = 5000;

	// Vinf at departure in ruv coordinates
	LB[1] = 0; // km/s
	UB[1] = 3; // km/s
	LB[2] = 0;
	UB[2] = 1;
	LB[3] = 0;
	UB[3] = 1;
	for (int i = 0; i < n; ++i) {
		// Segments thrust throttle
		LB[3 * i + 4] = 0;
		UB[3 * i + 4] = 1;
		// Segments thrust uv variables
		LB[3 * i + 5] = 0;
		UB[3 * i + 5] = 1;
		LB[3 * i + 6] = 0;
		UB[3 * i + 6] = 1;
	}
	// Transfer tim in days
	LB.back() = 0;
	UB.back() = 1000;
}

void earth_mars_lt2::human_readable(const std::vector<double> &x) const
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
	Dr = sqrt(Dr) * R;
	Dv = sqrt(Dv) * V;
	cout << endl << setw(25) << "Pos mismatch (km):" << right << setw(15) << Dr << left;// << setw(15) << (r_fwd[0] - r_back[0]) << setw(15) << (r_fwd[1] - r_back[1]) << setw(15) << (r_fwd[2] - r_back[2]) << left << endl;
	cout << endl << setw(25) << "Vel mismatch (km/sec):" << right << setw(15) << Dv << left << endl << endl;// << setw(15) << (v_fwd[0] - v_back[0]) << setw(15) << (v_fwd[1] - v_back[1]) << setw(15) << (v_fwd[2] - v_back[2]) << left << endl << endl;

	cout << setw(40) << "Time of flight:" << x.back() << " days" << endl;
	cout << setw(40) << "Total DV:" << main_objfun(x) * V << " Km/s" << endl;
	cout << setw(40) << "Final mass:" << M * std::exp(-main_objfun(x) * V * 1000. / (9.80665 * Isp)) << " Kg" << endl;
}

double earth_mars_lt2::objfun_(const std::vector<double> &x) const
{
	double r_fwd[3], v_fwd[3], r_back[3], v_back[3];
	state_mismatch(x,r_fwd,v_fwd,r_back,v_back);
	const double s_mismatch = std::sqrt((r_back[0] - r_fwd[0]) * (r_back[0] - r_fwd[0]) + (r_back[1] - r_fwd[1]) * (r_back[1] - r_fwd[1]) +
	                                    (r_back[2] - r_fwd[2]) * (r_back[2] - r_fwd[2]) + (v_back[0] - v_fwd[0]) * (v_back[0] - v_fwd[0]) + (v_back[1] - v_fwd[1]) * (v_back[1] - v_fwd[1]) +
	                                    (v_back[2] - v_fwd[2]) * (v_back[2] - v_fwd[2]));
	return main_objfun(x) + 1000 * s_mismatch;
}

double earth_mars_lt2::main_objfun(const std::vector<double> &x) const
{
	double retval = 0;
	double dt = (x.back() / n) * 86400 / T;
	for (int i = 0; i < n; ++i) {
		retval += x[3 * i + 4] * thrust / M * dt ;
	}
	return retval;
}

std::vector<std::vector<double> > earth_mars_lt2::visualize(const std::vector<double> &x, const int N) const
{
	std::vector<std::vector<double> > retval_fwd,retval_back;
	std::vector<double> xyz(10);

	const int n_seg_fwd = (n + 1) / 2, n_seg_back = n / 2;
	const double dt = (x.back() / n) * 86400 / T;

	double r_fwd[3], v_fwd[3], r_back[3], v_back[3];
	double fixed_thrust[3];

	//1 - Evaluate xyzdxdydzuxuyuz and store data
	earth_eph(x[0],r_fwd,v_fwd);
	kick(v_fwd,&x[1]);
	ruv2cart(fixed_thrust,&x[4]);

	//from throttle to non dimensional thrust
	fixed_thrust[0] *= thrust / M;
	fixed_thrust[1] *= thrust / M;
	fixed_thrust[2] *= thrust / M;
	//store
	for (int i = 0; i<3; ++i) {
		xyz[0] = 0;
		xyz[i+1] = r_fwd[i];
		xyz[i+4] = v_fwd[i];
		xyz[i+7] = fixed_thrust[i];
	}
	;
	retval_fwd.push_back(xyz);

	//Forward Propagation
	for (int i = 0; i < n_seg_fwd; ++i) {
		ruv2cart(fixed_thrust,&x[3*i+4]);   //We extract from the decision vector the throttle
		fixed_thrust[0] *= thrust / M;	    //And we evaluate the actual thrust (non dimensional)
		fixed_thrust[1] *= thrust / M;
		fixed_thrust[2] *= thrust / M;
		for (int j =0; j<N; ++j) {
			propagate(r_fwd,v_fwd,dt / (double)N, fixed_thrust);
			xyz[0] += dt / (double)N;
			for (int i = 0; i<3; ++i) {
				xyz[i+1] = r_fwd[i];
				xyz[i+4] = v_fwd[i];
				xyz[i+7] = fixed_thrust[i];
			};
			retval_fwd.push_back(xyz);
		}
	}

	mars_eph(x[0] + x.back(),r_back,v_back);
	ruv2cart(fixed_thrust,&x[x.size()-4]);
	fixed_thrust[0] *= thrust / M;
	fixed_thrust[1] *= thrust / M;
	fixed_thrust[2] *= thrust / M;

	xyz[0] = x.back() * 86400 / T;
	for (int i = 0; i<3; ++i) {
		xyz[i+1] = r_back[i];
		xyz[i+4] = v_back[i];
		xyz[i+7] = fixed_thrust[i];
	}
	;
	retval_back.push_back(xyz);

	//Bacward Propagation
	for (int i = 0; i < n_seg_back; ++i) {
		ruv2cart(fixed_thrust,&x[ (x.size()-1) - 3 * (i+1) ]);   //We extract from the decision vector the throttle
		fixed_thrust[0] *= thrust / M;	   //And we evaluate the actual thrust (non dimensional)
		fixed_thrust[1] *= thrust / M;
		fixed_thrust[2] *= thrust / M;
		for (int j =0; j<N; ++j) {
			back_propagate(r_back,v_back,dt / (double)N, fixed_thrust);
			xyz[0] -= dt / (double)N;
			for (int i = 0; i<3; ++i) {
				xyz[i+1] = r_back[i];
				xyz[i+4] = v_back[i];
				xyz[i+7] = fixed_thrust[i];
			}
			;
			retval_back.push_back(xyz);
		}
	}

	//Assemble the final result in the correct time progression
	std::reverse(retval_back.begin(), retval_back.end());
	for (int i=0; i<retval_back.size(); ++i) {
		retval_fwd.push_back(retval_back[i]);
	}

	return retval_fwd;
}

void earth_mars_lt2::state_mismatch(const std::vector<double> &x, double *r_fwd, double *v_fwd, double *r_back, double *v_back) const
{
	//We divide the n segments in forward segments and backward segments.
	const int n_seg_fwd = (n + 1) / 2, n_seg_back = n / 2;
	//This is the segment duration (in non dimensional units)
	const double dt = (x.back() / n) * 86400 / T;
	double fixed_thrust[3];
	//Here we evaluate the Earth position and velocity and copy these numbers into the forward spacecraft state
	earth_eph(x[0],r_fwd,v_fwd);
	//We add to the spacecraft velocity of the launcher DV
	kick(v_fwd,&x[1]);
	//Here we evaluate Mars position and velocity and copy these numbers into the backward spacecraft state
	mars_eph(x[0] + x.back(),r_back,v_back);

	//Forward Propagation
	for (int i = 0; i < n_seg_fwd; ++i) {
		ruv2cart(fixed_thrust,&x[3*i+4]);   //We extract from the decision vector the throttle
		fixed_thrust[0] *= thrust / M;	   //And we evaluate the actual thrust (non dimensional)
		fixed_thrust[1] *= thrust / M;
		fixed_thrust[2] *= thrust / M;
		propagate(r_fwd,v_fwd,dt,fixed_thrust);
	}

	//Backward Propagation
	for (int i = 0; i < n_seg_back; ++i) {
		ruv2cart(fixed_thrust,&x[ (x.size()-1) - 3 * (i+1) ]);
		fixed_thrust[0] *= thrust / M;
		fixed_thrust[1] *= thrust / M;
		fixed_thrust[2] *= thrust / M;
		back_propagate(r_back,v_back,dt,fixed_thrust);
	}
}

void earth_mars_lt2::ruv2cart(double *output, const double *input)
{
	double theta, phi, r = input[0];
	theta = 2 * M_PI * input[1];
	phi = std::acos(2 * input[2] - 1);
	output[0] = r * std::cos(theta) * std::sin(phi);
	output[1] = r * std::sin(theta) * std::sin(phi);
	output[2] = r * std::cos(phi);
}

void earth_mars_lt2::earth_eph(const double &mjd2000, double *position, double *velocity)
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

void earth_mars_lt2::mars_eph(const double &mjd2000, double *position, double *velocity)
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

void earth_mars_lt2::kick(double *output, const double *input)
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

void earth_mars_lt2::propagate(double *r, double *v, const double &t, double* fixed_thrust) const
{
	double y[6];

	//Set initial conditions
	y[0] = r[0];
	y[1] = r[1];
	y[2] = r[2];
	y[3] = v[0];
	y[4] = v[1];
	y[5] = v[2];

	//Integrate
	taylor_fixedthrust(y,0,t,fixed_thrust,eps,eps);

	//Copy state y into output r,v
	r[0] = y[0];
	r[1] = y[1];
	r[2] = y[2];
	v[0] = y[3];
	v[1] = y[4];
	v[2] = y[5];
}

void earth_mars_lt2::back_propagate(double *r, double *v, const double &t, double* fixed_thrust) const
{
	double y[6];

	//Set initial conditions
	y[0] = r[0];
	y[1] = r[1];
	y[2] = r[2];
	y[3] = -v[0];
	y[4] = -v[1];
	y[5] = -v[2];

	//Integrate
	taylor_fixedthrust(y,0,t,fixed_thrust,eps,eps);

	//Copy state y into output r,v
	r[0] = y[0];
	r[1] = y[1];
	r[2] = y[2];
	v[0] = -y[3];
	v[1] = -y[4];
	v[2] = -y[5];
}

}
}
