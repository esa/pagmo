/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#include <cmath>
#include "Astro_Functions.h"
#include "Lambert.h"
#include "mga_dsm.h"
#include "propagateKEP.h"
#include "time2distance.h"

using namespace std;

// [MR] TODO: exctract these somewhere...
const double MU[9] = {
	1.32712428e11,          // SUN                                  = 0
	22321,					// Gravitational constant of Mercury	= 1
	324860,					// Gravitational constant of Venus		= 2
	398601.19,				// Gravitational constant of Earth		= 3
	42828.3,				// Gravitational constant of Mars		= 4
	126.7e6,				// Gravitational constant of Jupiter	= 5
	0.37939519708830e8,		// Gravitational constant of Saturn		= 6
	5.78e6,					// Gravitational constant of Uranus		= 7
	6.8e6					// Gravitational constant of Neptune	= 8
};

//  Definition of planetari radii
//  TODO: maybe put missing values here so that indices correspond to those from MU[]?
const double RPL[6] = {
	2440,   // Mercury
	6052,   // Venus
	6378,   // Earth
	3397,   // Mars
	71492,  // Jupiter
	60330   // Saturn
};

//TODO: move somewhere else
void vector_normalize(const double in[3], double out[3])
{
	double norm = norm2(in);
	for(int i = 0; i < 3; i++) {
		out[i] = in[i] / norm;
	}
}

/**
 * Compute velocity and position of an celestial object of interest at specified time.
 *
 * problem - concerned problem
 * T       - time
 * i_count - hop number (starting from 0)
 * r       - [output] object's position
 * v       - [output] object's velocity
 */
void get_celobj_r_and_v(const mgadsmproblem& problem, const double T, const int i_count, double* r, double* v)
{
	if (problem.sequence[i_count] < 10) { //normal planet
		Planet_Ephemerides_Analytical (T, problem.sequence[i_count],
			r, v); // r and  v in heliocentric coordinate system
	} else { //asteroid
		Custom_Eph(T + 2451544.5, problem.asteroid.epoch, problem.asteroid.keplerian,
			r, v);
	}
}

/**
 * Precomputes all velocities and positions of celestial objects of interest for the problem.
 * Before calling this function, r and v verctors must be pre-allocated with sufficient amount of entries.
 *
 * problem - concerned problem
 * r       - [output] array of position vectors
 * v       - [output] array of velocity vectors
 */
void precalculate_ers_and_vees(const vector<double>& t, const mgadsmproblem& problem, std::vector<double*>& r, std::vector<double*>& v)
{
	double T = t[0]; //time of departure

	for(unsigned int i_count = 0; i_count < problem.sequence.size(); i_count++) {
		get_celobj_r_and_v(problem, T, i_count, r[i_count], v[i_count]);
		T += t[4 + i_count]; //time of flight
	}
}

/**
 * Get gravitational constant of an celestial object of interest.
 *
 * problem - concerned problem
 * i_count - hop number (starting from 0)
 */
double get_celobj_mu(const mgadsmproblem& problem, const int i_count)
{
	if (problem.sequence[i_count] < 10) { //normal planet
		return MU[problem.sequence[i_count]];
	} else { //asteroid
		return problem.asteroid.mu;
	}
}

// FIRST BLOCK (P1 to P2)
/**
 * t          - decision vector
 * problem    - problem parameters
 * r          - planet positions
 * v          - planet velocities
 * DV         - [output] velocity contributions table
 * v_sc_pl_in - [output] next hop input speed
 */
void first_block(const vector<double>& t, const mgadsmproblem& problem, const std::vector<double*>& r, std::vector<double*>& v, std::vector<double>& DV, double v_sc_nextpl_in[3])
{
	//First, some helper constants to make code more readable
	const int n = problem.sequence.size();
	const double VINF = t[1];         // Hyperbolic escape velocity (km/sec)
	const double udir = t[2];         // Hyperbolic escape velocity var1 (non dim)
	const double vdir = t[3];         // Hyperbolic escape velocity var2 (non dim)
	// [MR] {LITTLE HACKER TRICK} Instead of copying (!) arrays let's just introduce pointers to appropriate positions in the decision vector.
	const double *tof = &t[4];
	const double *alpha = &t[n+3];

	int i; //loop counter

	// Spacecraft position and velocity at departure
	double vtemp[3];
	cross(r[0], v[0], vtemp);

	double zP1[3];
	vector_normalize(vtemp, zP1);

	double iP1[3];
	vector_normalize(v[0], iP1);

	double jP1[3];
	cross(zP1, iP1, jP1);

	double theta, phi;
	theta = 2 * M_PI * udir;             // See Picking a Point on a Sphere
	phi = acos(2 * vdir - 1) - M_PI / 2; // In this way: -pi/2<phi<pi/2 so phi can be used as out-of-plane rotation

	double vinf[3];
	for (i = 0; i < 3; i++)
		vinf[i] = VINF * (cos(theta) * cos(phi) * iP1[i] + sin(theta) * cos(phi) * jP1[i] + sin(phi) * zP1[i]);

	//double v_sc_pl_in[3];  // Spacecraft absolute incoming velocity at P1 [MR] not needed?
	double v_sc_pl_out[3]; // Spacecraft absolute outgoing velocity at P1

	for (i = 0; i < 3; i++)
	{
		//v_sc_pl_in[i]  = v[0][i];
		v_sc_pl_out[i] = v[0][i] + vinf[i];
	}

	// Computing S/C position and absolute incoming velocity at DSM1
	double rd[3], v_sc_dsm_in[3];

	propagateKEP(r[0], v_sc_pl_out, alpha[0] * tof[0] * 86400, MU[0],
			rd, v_sc_dsm_in); // [MR] last two are output.

	// Evaluating the Lambert arc from DSM1 to P2
	double Dum_Vec[3]; // [MR] Rename it to something sensible...
	vett(rd, r[1], Dum_Vec);

	int lw = (Dum_Vec[2] > 0) ? 0 : 1;
	double a, p, theta2;
	int iter_unused; // [MR] unused variable

	double v_sc_dsm_out[3]; // DSM output speed

	LambertI(rd, r[1], tof[0] * (1 - alpha[0]) * 86400, MU[0], lw,
		v_sc_dsm_out, v_sc_nextpl_in, a, p, theta2, iter_unused);	// [MR] last 6 are output

	// First Contribution to DV (the 1st deep space maneuver)
	for (i = 0; i < 3; i++)
	{
		Dum_Vec[i] = v_sc_dsm_out[i] - v_sc_dsm_in[i]; // [MR] Temporary variable reused. Dirty.
	}

	DV[0] = norm2(Dum_Vec);
}

// ------
// INTERMEDIATE BLOCK
// WARNING: i_count starts from 0
double intermediate_block(const vector<double>& t, const mgadsmproblem& problem, const std::vector<double*>& r, const std::vector<double*>& v, int i_count, const double v_sc_pl_in[], std::vector<double>& DV, double* v_sc_nextpl_in)
{
	//[MR] A bunch of helper variables to simplify the code
	const int n = problem.sequence.size();
	// [MR] {LITTLE HACKER TRICK} Instead of copying (!) arrays let's just introduce pointers to appropriate positions in the decision vector.
	const double *tof = &t[4];
	const double *alpha = &t[n+3];
	const double *rp_non_dim = &t[2*n+2]; // non-dim perigee fly-by radius of planets P2..Pn(-1) (i=1 refers to the second planet)
	const double *gamma = &t[3*n];        // rotation of the bplane-component of the swingby outgoing
	const vector<int>& sequence = problem.sequence;

	int i; //loop counter

	// Evaluation of the state immediately after Pi
	double v_rel_in[3];
	double vrelin = 0.0;

	for (i = 0; i < 3; i++)
	{
		v_rel_in[i] = v_sc_pl_in[i] - v[i_count+1][i];
		vrelin += v_rel_in[i] * v_rel_in[i];
	}

	// Hop object's gravitional constant
	double hopobj_mu = get_celobj_mu(problem, i_count + 1);

	double e = 1.0 + rp_non_dim[i_count] * RPL[sequence[i_count + 1] - 1] * vrelin / hopobj_mu;

	double beta_rot = 2 * asin(1 / e); // velocity rotation

	double ix[3];
	vector_normalize(v_rel_in, ix);

	double vpervnorm[3];
	vector_normalize(v[i_count+1], vpervnorm);

	double iy[3];
	vett(ix, vpervnorm, iy);
	vector_normalize(iy, iy); // [MR]this *might* not work properly...

	double iz[3];
	vett(ix, iy, iz);

	double v_rel_in_norm = norm2(v_rel_in);

	double v_sc_pl_out[3]; // TODO: document me!

	for (i = 0; i < 3; i++)
	{
		double iVout = cos(beta_rot) * ix[i] + cos(gamma[i_count]) * sin(beta_rot) * iy[i] + sin(gamma[i_count]) * sin(beta_rot) * iz[i];
		double v_rel_out = v_rel_in_norm * iVout;
		v_sc_pl_out[i] = v[i_count + 1][i] + v_rel_out;
	}

	// Computing S/C position and absolute incoming velocity at DSMi
	double rd[3], v_sc_dsm_in[3];

	propagateKEP(r[i_count + 1], v_sc_pl_out, alpha[i_count+1] * tof[i_count+1] * 86400, MU[0],
					rd, v_sc_dsm_in); // [MR] last two are output

	// Evaluating the Lambert arc from DSMi to Pi+1
	double Dum_Vec[3]; // [MR] Rename it to something sensible...
	vett(rd, r[i_count + 2], Dum_Vec);

	int lw = (Dum_Vec[2] > 0) ? 0 : 1;
	double a, p, theta;
	int iter_unused; // [MR] unused variable

	double v_sc_dsm_out[3]; // DSM output speed

	LambertI(rd, r[i_count + 2], tof[i_count + 1] * (1 - alpha[i_count + 1]) * 86400, MU[0], lw,
		v_sc_dsm_out, v_sc_nextpl_in, a, p, theta, iter_unused); // [MR] last 6 are output.

	// DV contribution
	for (i = 0; i < 3; i++) {
		Dum_Vec[i] = v_sc_dsm_out[i] - v_sc_dsm_in[i]; // [MR] Temporary variable reused. Dirty.
	}

	DV[i_count + 1] = norm2(Dum_Vec);

	return vrelin;
}

// FINAL BLOCK
//
void final_block(const mgadsmproblem& problem, const std::vector<double*>& , const std::vector<double*>& v, const double v_sc_pl_in[], std::vector<double>& DV)
{
	//[MR] A bunch of helper variables to simplify the code
	const int n = problem.sequence.size();
	const double rp_target =  problem.rp;
	const double e_target = problem.e;
	const vector<int>& sequence = problem.sequence;

	int i; //loop counter

	// Evaluation of the arrival DV
	double Dum_Vec[3];
	for (i = 0; i < 3; i++) {
		Dum_Vec[i] = v[n-1][i] - v_sc_pl_in[i];
	}

	double DVrel, DVarr;
	DVrel = norm2(Dum_Vec);  // Relative velocity at target planet

	if ((problem.type  == orbit_insertion) || (problem.type == total_DV_orbit_insertion)) {
		double DVper = sqrt(DVrel * DVrel + 2 * MU[sequence[n - 1]] / rp_target); //[MR] should MU be changed to get_... ?
	    double DVper2 = sqrt(2 * MU[sequence[n - 1]] / rp_target - MU[sequence[n - 1]] / rp_target * (1 - e_target));
	    DVarr = fabs(DVper - DVper2);
	} else if (problem.type == rndv){
		DVarr = DVrel;
	} else if (problem.type == total_DV_rndv){
		DVarr = DVrel;
	} else {
		DVarr = 0.0;  // no DVarr is considered
	}

	DV[n - 1] = DVarr;
}


int MGA_DSM(
			/* INPUT values: */ //[MR] make this parameters const, if they are not modified and possibly references (especially 'problem').
			const vector<double> &t,	// it is the vector which provides time in modified julian date 2000. [MR] ??? Isn't it the decision vetor ???
			const mgadsmproblem& problem,

			/* OUTPUT values: */
			double &J    // output
			)
{
	//[MR] A bunch of helper variables to simplify the code
	const int n = problem.sequence.size();

	int i; //loop counter

	//References to objects pre-allocated in the mgadsm struct
	std::vector<double*>& r = problem.r;
	std::vector<double*>& v = problem.v;

	std::vector<double>& DV = problem.DV; //DV contributions

	precalculate_ers_and_vees(t, problem, r, v);

	double inter_pl_in_v[3], inter_pl_out_v[3]; //inter-hop velocities

	// FIRST BLOCK
	first_block(t, problem, r, v,
		DV, inter_pl_out_v); // [MR] output

	// INTERMEDIATE BLOCK

	for (int i_count=0; i_count < n - 2; i_count++)	{
		//copy previous output velocity to current input velocity
		inter_pl_in_v[0] = inter_pl_out_v[0]; inter_pl_in_v[1] = inter_pl_out_v[1]; inter_pl_in_v[2] = inter_pl_out_v[2];

		problem.vrelin_vec[i_count] = intermediate_block(t, problem, r, v, i_count, inter_pl_in_v,DV, inter_pl_out_v);
	}

	//copy previous output velocity to current input velocity
	inter_pl_in_v[0] = inter_pl_out_v[0]; inter_pl_in_v[1] = inter_pl_out_v[1]; inter_pl_in_v[2] = inter_pl_out_v[2];
	// FINAL BLOCK
	final_block(problem, r, v, inter_pl_in_v,
		DV);

	// **************************************************************************
	// Evaluation of total DV spent by the propulsion system
	// **************************************************************************
	double DVtot = 0.0;

	for (i = 0; i < n; i++) {
		DVtot += DV[i];
	}

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//[MR] Calculation of the actual procedure output (DVvec and J)

	const double& VINF = t[1];         // Variable renaming: Hyperbolic escape velocity (km/sec)

	for (i = n; i > 0; i--) {
		DV[i] = DV[i - 1];
	}
	DV[0] = VINF;


	// Finally our objective function (J) is:

	if (problem.type == orbit_insertion) {
		J = DVtot;
	} else if (problem.type == total_DV_orbit_insertion) {
		J = DVtot + VINF;
	} else if (problem.type == rndv) {
		J = DVtot;
	} else if (problem.type == total_DV_rndv) {
		J = DVtot + VINF;
	} else if (problem.type == time2AUs) { // [MR] TODO: extract method
		// [MR] helper constants
		const vector<int>& sequence = problem.sequence;
		const double *rp_non_dim = &t[2*n+2]; // non-dim perigee fly-by radius of planets P2..Pn(-1) (i=1 refers to the second planet)
		const double *gamma = &t[3*n];        // rotation of the bplane-component of the swingby outgoing
		const double AUdist = problem.AUdist;
		const double DVtotal = problem.DVtotal;
		const double DVonboard = problem.DVonboard;
		const double *tof = &t[4];

		// non dimensional units
		const double AU  = 149597870.66;
		const double V = sqrt(MU[0] / AU);
		const double T = AU / V;

		//evaluate the state of the spacecraft after the last fly-by
		double vrelin = 0.0;
		double v_rel_in[3];
		for (i=0; i<3; i++)
		{
			v_rel_in[i] = inter_pl_in_v[i] - v[n-1][i];
			vrelin += v_rel_in[i] * v_rel_in[i];
		}

		double e = 1.0 + rp_non_dim[n - 2] * RPL[sequence[n - 2 + 1] - 1] * vrelin / get_celobj_mu(problem, n - 1); //I hope the planet index (n - 1) is OK

		double beta_rot=2*asin(1/e);              // velocity rotation

		double vrelinn = norm2(v_rel_in);
		double ix[3];
		for (i=0; i<3; i++)
			ix[i] = v_rel_in[i]/vrelinn;
			// ix=r_rel_in/norm(v_rel_in);  // activating this line and disactivating the one above
                                 // shifts the singularity for r_rel_in parallel to v_rel_in

		double vnorm = norm2(v[n-1]);

        double vtemp[3];
		for (i=0; i<3; i++)
			vtemp[i] = v[n-1][i]/vnorm;

		double iy[3];
		vett(ix, vtemp, iy);

		double iynorm = norm2(iy);
		for (i=0; i<3; i++)
			iy[i] = iy[i]/iynorm;

		double iz[3];
		vett(ix, iy, iz);
		double v_rel_in_norm = norm2(v_rel_in);

		double v_sc_pl_out[3]; // TODO: document me!
		for (i = 0; i < 3; i++)
		{
			double iVout = cos(beta_rot) * ix[i] + cos(gamma[n - 2]) * sin(beta_rot) * iy[i] + sin(gamma[n - 2]) * sin(beta_rot) * iz[i];
			double v_rel_out = v_rel_in_norm * iVout;
			v_sc_pl_out[i] = v[n - 1][i] + v_rel_out;
		}

		double r_per_AU[3];
		double v_sc_pl_out_per_V[3];
		for (i = 0; i < 3; i++)
		{
			r_per_AU[i] = r[n - 1][i] / AU;
			v_sc_pl_out_per_V[i] = v_sc_pl_out[i] / V;
		}

		double time = time2distance(r_per_AU, v_sc_pl_out_per_V, AUdist);
		// if (time == -1) cout << " ERROR" << endl;

		if (time != -1)
		{
			double DVpen=0;
			double sum = 0.0;

			for (i=0; i<n+1; i++)
				sum += DV[i];

			if (sum > DVtotal)
				DVpen += DVpen+(sum-DVtotal);

			sum = 0.0;
			for (i=1; i<n+1; i++)
				sum+=DV[i];

			if (sum > DVonboard)
				DVpen = DVpen + (sum - DVonboard);

			sum = 0.0;
			for (i=0; i<n-1; i++)
				sum += tof[i];

			J= (time*T/60/60/24 + sum)/365.25 + DVpen*100;
		}
		else
			J = 100000;   // there was an ERROR in time2distance
	} // time2AU

	return 0;
}
