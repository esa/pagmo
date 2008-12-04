// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //

#include <iostream>
#include <iomanip>
#include "Pl_Eph_An.h"
#include "Astro_Functions.h"

void Planet_Ephemerides_Analytical (const double mjd2000,
									const int planet,
									double *position,
									double *velocity)
{
	const double pi = acos(-1.0);
	const double RAD = pi/180.0;
	const double AU = 149597870.66; // Astronomical Unit
	const double KM = AU;
 	//const double MuSun = 1.32712440018e+11; //Gravitational constant of Sun);
	const double MuSun = 1.327124280000000e+011;  //Gravitational constant of Sun);
	double Kepl_Par[6];
	double XM;

	double T =(mjd2000 + 36525.00)/36525.00;

	switch (planet)
	{
		case(1):// Mercury
			Kepl_Par[0]=(0.38709860);
			Kepl_Par[1]=(0.205614210 + 0.000020460*T - 0.000000030*T*T);
			Kepl_Par[2]=(7.002880555555555560 + 1.86083333333333333e-3*T - 1.83333333333333333e-5*T*T);
			Kepl_Par[3]=(4.71459444444444444e+1 + 1.185208333333333330*T + 1.73888888888888889e-4*T*T);
			Kepl_Par[4]=(2.87537527777777778e+1 + 3.70280555555555556e-1*T +1.20833333333333333e-4*T*T);
			XM   = 1.49472515288888889e+5 + 6.38888888888888889e-6*T;
			Kepl_Par[5]=(1.02279380555555556e2 + XM*T);
		break;
		case(2):// Venus
			Kepl_Par[0]=(0.72333160);
			Kepl_Par[1]=(0.006820690 - 0.000047740*T + 0.0000000910*T*T);
			Kepl_Par[2]=(3.393630555555555560 + 1.00583333333333333e-3*T - 9.72222222222222222e-7*T*T);
			Kepl_Par[3]=(7.57796472222222222e+1 + 8.9985e-1*T + 4.1e-4*T*T);
			Kepl_Par[4]=(5.43841861111111111e+1 + 5.08186111111111111e-1*T -1.38638888888888889e-3*T*T);
			XM =5.8517803875e+4 + 1.28605555555555556e-3*T;
			Kepl_Par[5]=(2.12603219444444444e2 + XM*T);
		break;
		case(3):// Earth
			Kepl_Par[0]=(1.000000230);
			Kepl_Par[1]=(0.016751040 - 0.000041800*T - 0.0000001260*T*T);
			Kepl_Par[2]=(0.00);
			Kepl_Par[3]=(0.00);
			Kepl_Par[4]=(1.01220833333333333e+2 + 1.7191750*T + 4.52777777777777778e-4*T*T + 3.33333333333333333e-6*T*T*T);
			 XM   = 3.599904975e+4 - 1.50277777777777778e-4*T - 3.33333333333333333e-6*T*T;
			Kepl_Par[5]=(3.58475844444444444e2 + XM*T);
		break;
		case(4):// Mars
			Kepl_Par[0]=(1.5236883990);
			Kepl_Par[1]=(0.093312900 + 0.0000920640*T - 0.0000000770*T*T);
			Kepl_Par[2]=(1.850333333333333330 - 6.75e-4*T + 1.26111111111111111e-5*T*T);
			Kepl_Par[3]=(4.87864416666666667e+1 + 7.70991666666666667e-1*T - 1.38888888888888889e-6*T*T - 5.33333333333333333e-6*T*T*T);
			Kepl_Par[4]=(2.85431761111111111e+2 + 1.069766666666666670*T +  1.3125e-4*T*T + 4.13888888888888889e-6*T*T*T);
			XM   = 1.91398585e+4 + 1.80805555555555556e-4*T + 1.19444444444444444e-6*T*T;
			Kepl_Par[5]=(3.19529425e2 + XM*T);
		break;
		case(5):// Jupiter
			Kepl_Par[0]=(5.2025610);
			Kepl_Par[1]=(0.048334750 + 0.000164180*T  - 0.00000046760*T*T -0.00000000170*T*T*T);
			Kepl_Par[2]=(1.308736111111111110 - 5.69611111111111111e-3*T +  3.88888888888888889e-6*T*T);
			Kepl_Par[3]=(9.94433861111111111e+1 + 1.010530*T + 3.52222222222222222e-4*T*T - 8.51111111111111111e-6*T*T*T);
			Kepl_Par[4]=(2.73277541666666667e+2 + 5.99431666666666667e-1*T + 7.0405e-4*T*T + 5.07777777777777778e-6*T*T*T);
			XM   = 3.03469202388888889e+3 - 7.21588888888888889e-4*T + 1.78444444444444444e-6*T*T;
			Kepl_Par[5]=(2.25328327777777778e2 + XM*T);
		break;
		case(6):// Saturn
			Kepl_Par[0]=(9.5547470);
			Kepl_Par[1]=(0.055892320 - 0.00034550*T - 0.0000007280*T*T + 0.000000000740*T*T*T);
			Kepl_Par[2]=(2.492519444444444440 - 3.91888888888888889e-3*T - 1.54888888888888889e-5*T*T + 4.44444444444444444e-8*T*T*T);
			Kepl_Par[3]=(1.12790388888888889e+2 + 8.73195138888888889e-1*T -1.52180555555555556e-4*T*T - 5.30555555555555556e-6*T*T*T);
			Kepl_Par[4]=(3.38307772222222222e+2 + 1.085220694444444440*T + 9.78541666666666667e-4*T*T + 9.91666666666666667e-6*T*T*T);
			XM   = 1.22155146777777778e+3 - 5.01819444444444444e-4*T - 5.19444444444444444e-6*T*T;
			Kepl_Par[5]=(1.75466216666666667e2 + XM*T);
		break;
		case(7):// Uranus
			Kepl_Par[0]=(19.218140);
			Kepl_Par[1]=(0.04634440 - 0.000026580*T + 0.0000000770*T*T);
			Kepl_Par[2]=(7.72463888888888889e-1 + 6.25277777777777778e-4*T + 3.95e-5*T*T);
			Kepl_Par[3]=(7.34770972222222222e+1 + 4.98667777777777778e-1*T + 1.31166666666666667e-3*T*T);
			Kepl_Par[4]=(9.80715527777777778e+1 + 9.85765e-1*T - 1.07447222222222222e-3*T*T - 6.05555555555555556e-7*T*T*T);
			XM   = 4.28379113055555556e+2 + 7.88444444444444444e-5*T + 1.11111111111111111e-9*T*T;
			Kepl_Par[5]=(7.26488194444444444e1 + XM*T);
		break;
		case(8)://Neptune
			Kepl_Par[0]=(30.109570);
			Kepl_Par[1]=(0.008997040 + 0.0000063300*T - 0.0000000020*T*T);
			Kepl_Par[2]=(1.779241666666666670 - 9.54361111111111111e-3*T - 9.11111111111111111e-6*T*T);
			Kepl_Par[3]=(1.30681358333333333e+2 + 1.0989350*T + 2.49866666666666667e-4*T*T - 4.71777777777777778e-6*T*T*T);
			Kepl_Par[4]=(2.76045966666666667e+2 + 3.25639444444444444e-1*T + 1.4095e-4*T*T + 4.11333333333333333e-6*T*T*T);
			XM   = 2.18461339722222222e+2 - 7.03333333333333333e-5*T;
			Kepl_Par[5]=(3.77306694444444444e1 + XM*T);
		break;
		case(9):// Pluto
			//Fifth order polynomial least square fit generated by Dario Izzo
			//(ESA ACT). JPL405 ephemerides (Charon-Pluto barycenter) have been used to produce the coefficients.
			//This approximation should not be used outside the range 2000-2100;
			T =mjd2000/36525.00;
			Kepl_Par[0]=(39.34041961252520 + 4.33305138120726*T - 22.93749932403733*T*T + 48.76336720791873*T*T*T - 45.52494862462379*T*T*T*T + 15.55134951783384*T*T*T*T*T);
			Kepl_Par[1]=(0.24617365396517 + 0.09198001742190*T - 0.57262288991447*T*T + 1.39163022881098*T*T*T - 1.46948451587683*T*T*T*T + 0.56164158721620*T*T*T*T*T);
			Kepl_Par[2]=(17.16690003784702 - 0.49770248790479*T + 2.73751901890829*T*T - 6.26973695197547*T*T*T + 6.36276927397430*T*T*T*T - 2.37006911673031*T*T*T*T*T);
			Kepl_Par[3]=(110.222019291707 + 1.551579150048*T - 9.701771291171*T*T + 25.730756810615*T*T*T - 30.140401383522*T*T*T*T + 12.796598193159 * T*T*T*T*T);
			Kepl_Par[4]=(113.368933916592 + 9.436835192183*T - 35.762300003726*T*T + 48.966118351549*T*T*T - 19.384576636609*T*T*T*T - 3.362714022614 * T*T*T*T*T);
			Kepl_Par[5]=(15.17008631634665 + 137.023166578486*T + 28.362805871736*T*T - 29.677368415909*T*T*T - 3.585159909117*T*T*T*T + 13.406844652829 * T*T*T*T*T);
		break;

	}

	// conversion of AU into KM
	Kepl_Par[0] *= KM;

	// conversion of DEG into RAD
	Kepl_Par[2] *= RAD;
	Kepl_Par[3] *= RAD;
	Kepl_Par[4] *= RAD;
	Kepl_Par[5] *= RAD;
	Kepl_Par[5] =  fmod(Kepl_Par[5], 2.0*pi);

	// Conversion from Mean Anomaly to Eccentric Anomaly via Kepler's equation
	Kepl_Par[5] = Mean2Eccentric(Kepl_Par[5], Kepl_Par[1]);

	// Position and Velocity evaluation according to j2000 system
	Conversion(Kepl_Par, position, velocity, MuSun);
}

void Custom_Eph(const double jd,
				const double epoch,
				const double keplerian[],
				double *position,
				double *velocity)
{
	const double pi = acos(-1.0);
	const double RAD = pi/180.0;
	const double AU = 149597870.66; // Astronomical Unit
	const double muSUN = 1.32712428e+11;    // Gravitational constant of Sun
	double a,e,i,W,w,M,jdepoch,DT,n,E;
	double V[6];

	a=keplerian[0]*AU; // in km
    e=keplerian[1];
    i=keplerian[2];
	W=keplerian[3];
    w=keplerian[4];
    M=keplerian[5];
	jdepoch = epoch +2400000.5;
    DT=(jd-jdepoch)*86400;
    n=sqrt(muSUN/pow(a,3));

    M=M/180.0*pi;
    M+=n*DT;
    M=fmod(M,2*pi);
    E=Mean2Eccentric(M,e);
	V[0] = a; V[1] = e; V[2] = i*RAD;
	V[3] = W*RAD; V[4] = w*RAD; V[5] = E;

	Conversion(V,position,velocity,muSUN);
}
