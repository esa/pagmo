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

#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <iostream>
#include "Pl_Eph_An.h"
#include "mga.h"
#include "Lambert.h"
#include "PowSwingByInv.h"
#include "Astro_Functions.h"
#define MAX(a, b) (a > b ? a : b)

using namespace std;

//the function return 0 if the input is right or -1 it there is something wrong

int MGA(vector<double> t,	// it is the vector which provides time in modified julian date 2000.
								// The first entry is launch date, the next entries represent the time needed to
								// fly from last swing-by to current swing-by.
			mgaproblem problem,

			/* OUTPUT values: */
			vector<double>& rp,  // periplanets radius
			vector<double>& DV,  // final delta-Vs
			double &obj_funct)   //objective function

{
	const int n = problem.sequence.size();
	const vector<int> sequence = problem.sequence;
	const vector<int> rev_flag = problem.rev_flag;// array containing 0 clockwise, 1 un-clockwise
	customobject cust_obj = problem.asteroid;

	double MU[9] = {//1.32712440018e11, //SUN = 0
					1.32712428e11,
					22321,		// Gravitational constant of Mercury	= 1
					324860,		// Gravitational constant of Venus		= 2
					398601.19,	// Gravitational constant of Earth		= 3
					42828.3,	// Gravitational constant of Mars		= 4
					126.7e6,	// Gravitational constant of Jupiter	= 5
					37.9e6,		// Gravitational constant of Saturn		= 6
					5.78e6,		// Gravitational constant of Uranus		= 7
					6.8e6		// Gravitational constant of Neptune	= 8
				    };
	double penalty[9] = {0,
		                0,        // Mercury
						6351.8,   // Venus
						6778.1,   // Earth
						6000,     // Mars
						//671492, // Jupiter
                        600000,   // Jupiter
						70000,    // Saturn
                        0,        // Uranus
						0         // Neptune
	};

	double penalty_coeffs[9] = {0,
								0,      // Mercury
								0.01,   // Venus
								0.01,   // Earth
								0.01,   // Mars
								0.001,  // Jupiter
								0.01,   // Saturn
								0,      // Uranus
								0       // Neptune
	};

	double DVtot = 0;
	double Dum_Vec[3],Vin,Vout;
	double V_Lamb[2][2][3],dot_prod;
	double a,p,theta,alfa;
	double DVrel, DVarr=0;

	//only used for orbit insertion (ex: cassini)
		double DVper, DVper2;
		const double rp_target = problem.rp;
		const double e_target = problem.e;
		const double DVlaunch = problem.DVlaunch;

	//only used for asteroid impact (ex: gtoc1)
		const double initial_mass = problem.mass;   // Satellite initial mass [Kg]
		double final_mass;							// satelite final mass
		const double Isp = problem.Isp;             // Satellite specific impulse [s]
		const double g = 9.80665 / 1000.0;          // Gravity




	double *vec, *rec;
	vector<double*> r;    // {0...n-1} position
	vector<double*> v;    // {0...n-1} velocity

	double T = 0.0;         // total time

	int i_count, j_count, lw;

	int iter = 0;

	if (n >= 2)
	{
		for ( i_count = 0; i_count < n; i_count++)
		{
			vec = new double [3];  // velocity and position are 3 D vector
			rec = new double [3];
			r.push_back(vec);
			v.push_back(rec);

			DV [i_count] = 0.0;
		}

		T = 0;
		for (i_count = 0; i_count < n; i_count++)
		{
			T += t[i_count];
			if (sequence[i_count]<10)
				Planet_Ephemerides_Analytical (T, sequence[i_count],
					r[i_count], v[i_count]); //r and  v in heliocentric coordinate system
			else
			{
				Custom_Eph(T+2451544.5, cust_obj.epoch, cust_obj.keplerian, r[i_count], v[i_count]);
			}
		}

		vett(r[0], r[1], Dum_Vec);

		if (Dum_Vec[2] > 0)
			lw = (rev_flag[0] == 0) ? 0 : 1;
	    else
			lw = (rev_flag[0] == 0) ? 1 : 0;

		LambertI(r[0],r[1],t[1]*24*60*60,MU[0],lw,          // INPUT
			     V_Lamb[0][0],V_Lamb[0][1],a,p,theta,iter); // OUTPUT
		DV[0] = norm(V_Lamb[0][0], v[0]);                // Earth launch

		for (i_count = 1; i_count <= n-2; i_count++)
		{
			vett(r[i_count], r[i_count+1], Dum_Vec);

			if (Dum_Vec[2] > 0)
				lw = (rev_flag[i_count] == 0) ? 0 : 1;
			else
				lw = (rev_flag[i_count] == 0) ? 1 : 0;

			/*if (i_count%2 != 0)	{*/
			LambertI(r[i_count],r[i_count+1],t[i_count + 1]*24*60*60,MU[0],lw, // INPUT
				   V_Lamb[1][0],V_Lamb[1][1],a,p,theta,iter);                  // OUTPUT

			// norm first perform the subtraction of vet1-vet2 and the evaluate ||...||
			Vin  = norm(V_Lamb[0][1], v[i_count]);
			Vout = norm(V_Lamb[1][0], v[i_count]);

			dot_prod = 0.0;
			for (int i = 0; i < 3; i++)
			{
				dot_prod += (V_Lamb[0][1][i] - v[i_count][i]) * (V_Lamb[1][0][i] - v[i_count][i]);
			}
			alfa = acos ( dot_prod /(Vin * Vout) );

			// calculation of delta V at pericenter
			PowSwingByInv(Vin, Vout, alfa, DV[i_count], rp[i_count - 1]);

			rp[i_count - 1] *= MU[sequence[i_count]];

			if (i_count != n-2)  //swap
		    	for (j_count = 0; j_count < 3; j_count++)
				{
					V_Lamb[0][0][j_count] = V_Lamb[1][0][j_count];  // [j_count];
					V_Lamb[0][1][j_count] = V_Lamb[1][1][j_count];  // [j_count];
				}
			}
	}
	else
	{
		return -1;
	}

	for (i_count = 0; i_count < 3; i_count++)
		Dum_Vec[i_count] = v[n-1][i_count] - V_Lamb[1][1][i_count];

	DVrel = norm2(Dum_Vec);

	if (problem.type == total_DV_orbit_insertion){

		DVper  = sqrt(DVrel*DVrel + 2*MU[sequence[n-1]]/rp_target);
		DVper2 = sqrt(2*MU[sequence[n-1]]/rp_target - MU[sequence[n-1]]/rp_target*(1-e_target));
		DVarr = fabs(DVper - DVper2);
	}

	else if (problem.type == asteroid_impact){

		DVarr = DVrel;
	}



	DVtot = 0;

	for (i_count = 1; i_count < n-1; i_count++)
		DVtot += DV[i_count];

	if (problem.type == total_DV_orbit_insertion){

		DVtot += DVarr;
	}



	// Build Penalty
	for (i_count = 0;i_count < n-2; i_count++)
		if (rp[i_count] < penalty[sequence[i_count+1]])
			DVtot += penalty_coeffs[sequence[i_count+1]]*fabs(rp[i_count] - penalty[sequence[i_count+1]]);

	// Launcher Constraint
	if (DV[0] > DVlaunch)
		DVtot += (DV[0] - DVlaunch);

	if (problem.type == total_DV_orbit_insertion){

		obj_funct = DVtot;
	}

	else if (problem.type == asteroid_impact){

		// Evaluation of satellite final mass
		obj_funct = final_mass = initial_mass * exp(- DVtot/ (Isp * g));

		// V asteroid - V satellite
		for (i_count = 0; i_count < 3; i_count++)
			Dum_Vec[i_count] = v[n-1][i_count] - V_Lamb[1][1][i_count];// arrival relative velocity at the asteroid;

		dot_prod = 0;
		for (i_count = 0; i_count < 3 ; i_count++)
			dot_prod += Dum_Vec[i_count] * v[n-1][i_count];

		obj_funct = - (final_mass)* fabs(dot_prod);
	}


	// final clean
	for ( i_count = 0;i_count < n;i_count++)
	{
		delete [] r[i_count];
		delete [] v[i_count];
	}
	r.clear();
	v.clear();
	return 0;
}



