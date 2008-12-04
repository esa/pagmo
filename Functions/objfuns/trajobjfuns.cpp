/*
 *  GOProblems.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/17/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include <math.h>
#include "mga.h"
#include "mga_dsm.h"
#include "misc4Tandem.h"


using namespace std;

double messenger(const vector<double>& x){
	double Delta_V[6];
	mgadsmproblem problem;

	int sequence_[5] = {3, 3, 2, 2, 1}; // sequence of planets
	problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_+ 5 );
	problem.type = total_DV_rndv;
	
	
	//Memory allocation
	problem.r = std::vector<double*>(5);
	problem.v = std::vector<double*>(5);
	problem.DV = std::vector<double>(5+1);
	for(int i = 0; i < 5; i++) {
		problem.r[i] = new double[3];
		problem.v[i] = new double[3];
	}
	

	double obj = 0;

	MGA_DSM(
			/* INPUT values: */
			x,
		    problem,

			/* OUTPUT values: */
			obj, Delta_V);
			
	
	//Memory release
	for(int i = 0; i < 5; i++) {
		delete[] problem.r[i];
		delete[] problem.v[i];
	}
	problem.r.clear();
	problem.v.clear();
	
	
	return obj;
}

double cassini2(const vector<double>& x){
	double Delta_V[7];
	mgadsmproblem problem;

	int sequence_[6] = {3, 2, 2, 3, 5, 6}; // sequence of planets
	problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_+ 6 );
	problem.type = total_DV_rndv;

	double obj = 0;
	

	//Allocate temporary memory for MGA_DSM
	problem.r = std::vector<double*>(6);
	problem.v = std::vector<double*>(6);
	problem.DV = std::vector<double>(6+1);
	
	for(int i = 0; i < 6; i++) {
		problem.r[i] = new double[3];
		problem.v[i] = new double[3];
	}
	

	MGA_DSM(
			/* INPUT values: */
			x,
		    problem,

			/* OUTPUT values: */
			obj, Delta_V);
			

	//Free temporary memory for MGA_DSM
	for(int i = 0; i < 6; i++) {
		delete[] problem.r[i];
		delete[] problem.v[i];
	}
	problem.r.clear();
	problem.v.clear();
	
	
	return obj;
}

double rosetta(const vector<double>& x){
	double Delta_V[7];
	mgadsmproblem problem;

	int sequence_[6] = {3, 3, 4, 3, 3, 10}; // sequence of planets
	problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_+ 6 );
	problem.type = rndv;
	problem.asteroid.keplerian[0] = 3.50294972836275;
	problem.asteroid.keplerian[1] = 0.6319356;
	problem.asteroid.keplerian[2] =  7.12723;
	problem.asteroid.keplerian[3] = 	50.92302;
	problem.asteroid.keplerian[4] =  11.36788;
	problem.asteroid.keplerian[5] = 0.0;
	problem.asteroid.epoch = 52504.23754000012;
	problem.asteroid.mu = 0.0;
	

	//Allocate temporary memory for MGA_DSM
	problem.r = std::vector<double*>(6);
	problem.v = std::vector<double*>(6);
	problem.DV = std::vector<double>(6+1);
	
	for(int i = 0; i < 6; i++) {
		problem.r[i] = new double[3];
		problem.v[i] = new double[3];
	}
	

	double obj = 0;

	MGA_DSM(
			/* INPUT values: */
			x,
		    problem,

			/* OUTPUT values: */
			obj, Delta_V);
			

	//Free temporary memory for MGA_DSM
	for(int i = 0; i < 6; i++) {
		delete[] problem.r[i];
		delete[] problem.v[i];
	}
	problem.r.clear();
	problem.v.clear();
			
			
	return obj;
}


double tandem(const vector<double>& x){
	const int seqlen = 5;
	const int sequence_[seqlen] = {3, 2, 3, 3, 6};		// sequence of planets
	double Delta_V[seqlen+1];
	double obj = 0;
	mgadsmproblem problem;

	//defining the problem
	problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_+ seqlen );
	problem.type = orbit_insertion;
	problem.rp = 80330.0;
	problem.e = 0.98531407996358;
	
	
	//Allocate temporary memory for MGA_DSM
	problem.r = std::vector<double*>(seqlen);
	problem.v = std::vector<double*>(seqlen);
	problem.DV = std::vector<double>(seqlen+1);
	
	for(int i = 0; i < seqlen; i++) {
		problem.r[i] = new double[3];
		problem.v[i] = new double[3];
	}


	//calling mgadsm
	MGA_DSM(x,problem,obj,Delta_V);

	//evaluating the mass from the dvs
	double rE[3];
	double vE[3];
	Planet_Ephemerides_Analytical (x[0],3,rE,vE);
	double VINFE = x[1];
	double udir = x[2];
	double vdir = x[3];
	double vtemp[3];
	vtemp[0]= rE[1]*vE[2]-rE[2]*vE[1];
	vtemp[1]= rE[2]*vE[0]-rE[0]*vE[2];
	vtemp[2]= rE[0]*vE[1]-rE[1]*vE[0];
	double iP1[3];
	double normvE=sqrt(vE[0]*vE[0]+vE[1]*vE[1]+vE[2]*vE[2]);
	iP1[0]=	vE[0]/normvE;
	iP1[1]=	vE[1]/normvE;
	iP1[2]=	vE[2]/normvE;
	double zP1[3];
	double normvtemp=sqrt(vtemp[0]*vtemp[0]+vtemp[1]*vtemp[1]+vtemp[2]*vtemp[2]);
	zP1[0]= vtemp[0]/normvtemp;
	zP1[1]= vtemp[1]/normvtemp;
	zP1[2]= vtemp[2]/normvtemp;
	double jP1[3];
	jP1[0]= zP1[1]*iP1[2]-zP1[2]*iP1[1];
	jP1[1]= zP1[2]*iP1[0]-zP1[0]*iP1[2];
	jP1[2]= zP1[0]*iP1[1]-zP1[1]*iP1[0];
	double theta=2*M_PI*udir; 		//See Picking a Point on a Sphere
	double phi=acos(2*vdir-1)-M_PI/2; //In this way: -pi/2<phi<pi/2 so phi can be used as out-of-plane rotation
	double vinf[3];
	vinf[0]=VINFE*(cos(theta)*cos(phi)*iP1[0]+sin(theta)*cos(phi)*jP1[0]+sin(phi)*zP1[0]);
	vinf[1]=VINFE*(cos(theta)*cos(phi)*iP1[1]+sin(theta)*cos(phi)*jP1[1]+sin(phi)*zP1[1]);
	vinf[2]=VINFE*(cos(theta)*cos(phi)*iP1[2]+sin(theta)*cos(phi)*jP1[2]+sin(phi)*zP1[2]);
	//We rotate it to the equatorial plane
	ecl2equ(vinf,vinf);
	//And we find the declination in degrees
	double normvinf=sqrt(vinf[0]*vinf[0]+vinf[1]*vinf[1]+vinf[2]*vinf[2]);
	double sindelta = vinf[2] / normvinf;
	double declination = asin(sindelta)/M_PI*180;

	//double m_initial = SoyuzFregat(VINFE,declination);
	double m_initial = Atlas501(VINFE,declination);

	//We evaluate the final mass
	double Isp = 312;
	double g0 = 9.80665;
	double sumDVvec=0;
	//double totaltime=x[4]+x[5]+x[6]+x[7];
	for(unsigned int i=1;i<=problem.sequence.size();i++) {
		sumDVvec=sumDVvec+Delta_V[i];
	}
	double m_final;
	sumDVvec=sumDVvec+0.165; //losses for 3 swgbys + insertion

	//if ((totaltime>3652.5)|| (totaltime<3287.25))  {
	//	m_final = 0;
	//}
	//else {
		m_final = m_initial * exp(-sumDVvec/Isp/g0*1000);
	//}
	//return (2000-m_final)/1000;
	//return -log(m_final/m_initial)*Isp*g0/1000;
	
	
	//Free temporary memory for MGA_DSM
	for(int i = 0; i < seqlen; i++) {
		delete[] problem.r[i];
		delete[] problem.v[i];
	}
	problem.r.clear();
	problem.v.clear();	
	
	
	return -log(m_final);
}

double gtoc1( const vector<double>& x)
{
	const int GTOC1_DIM = 8;
	vector<double> Delta_V(GTOC1_DIM);
	vector<double> rp(GTOC1_DIM-2);
	vector<double> t(GTOC1_DIM);
	mgaproblem problem;

	//Filling up the problem parameters
	problem.type = asteroid_impact;
	problem.mass = 1500.0;				// Satellite initial mass [Kg]
	problem.Isp = 2500.0;               // Satellite specific impulse [s]
	problem.DVlaunch = 2.5;				// Launcher DV in km/s

	int sequence_[GTOC1_DIM] = {3,2,3,2,3,5,6,10}; // sequence of planets
	vector<int> sequence(GTOC1_DIM);
	problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_+GTOC1_DIM);

	const int rev_[GTOC1_DIM] = {0,0,0,0,0,0,1,0}; // sequence of clockwise legs
	vector<int> rev(GTOC1_DIM);
	problem.rev_flag.insert(problem.rev_flag.begin(), rev_, rev_+GTOC1_DIM);

	problem.asteroid.keplerian[0] = 2.5897261;    // Asteroid data
	problem.asteroid.keplerian[1] = 0.2734625;
	problem.asteroid.keplerian[2] = 6.40734;
	problem.asteroid.keplerian[3] = 128.34711;
	problem.asteroid.keplerian[4] = 264.78691;
	problem.asteroid.keplerian[5] = 320.479555;
	problem.asteroid.epoch = 53600;

	double obj = 0;

	MGA(x,problem,rp,Delta_V,obj);

	return obj;
}

double cassini1( const vector<double>& x)
{
	const int CASSINI_DIM = 6;
	vector<double> Delta_V(CASSINI_DIM);
	vector<double> rp(CASSINI_DIM-2);
	vector<double> t(CASSINI_DIM);
	mgaproblem problem;

	//Filling up the problem parameters
	problem.type = total_DV_orbit_insertion;


	int sequence_[CASSINI_DIM] = {3,2,2,3,5,6}; // sequence of planets
	vector<int> sequence(CASSINI_DIM);
	problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_+CASSINI_DIM);

	const int rev_[CASSINI_DIM] = {0,0,0,0,0,0}; // sequence of clockwise legs
	vector<int> rev(CASSINI_DIM);
	problem.rev_flag.insert(problem.rev_flag.begin(), rev_, rev_+CASSINI_DIM);

	problem.e  = 0.98;      // Final orbit eccentricity
	problem.rp = 108950;    // Final orbit pericenter
	problem.DVlaunch = 0;   // Launcher DV

	double obj = 0;

	MGA(x,problem,rp,Delta_V,obj);

	return obj;
}

double sagas(const vector<double>& x){
	double Delta_V[4];
	mgadsmproblem problem;

	int sequence_[3] = {3,3,5}; // sequence of planets
	problem.sequence.insert(problem.sequence.begin(), sequence_, sequence_+ 3 );
	problem.type = time2AUs;
	problem.AUdist = 50.0;
	problem.DVtotal = 6.782;
	problem.DVonboard = 1.782;
	

	//Allocate temporary memory for MGA_DSM
	problem.r = std::vector<double*>(3);
	problem.v = std::vector<double*>(3);
	problem.DV = std::vector<double>(3+1);
	
	for(int i = 0; i < 3; i++) {
		problem.r[i] = new double[3];
		problem.v[i] = new double[3];
	}


	double obj = 0;

	MGA_DSM(
			/* INPUT values: */
			x,
		    problem,

			/* OUTPUT values: */
			obj, Delta_V);
			

	//Free temporary memory for MGA_DSM
	for(int i = 0; i < 3; i++) {
		delete[] problem.r[i];
		delete[] problem.v[i];
	}
	problem.r.clear();
	problem.v.clear();
			
	return obj;
}

