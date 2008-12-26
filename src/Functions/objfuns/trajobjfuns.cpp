/*
 *  GOProblems.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/17/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include <cmath>

#include "mga.h"
#include "mga_dsm.h"
#include "misc4Tandem.h"

using namespace std;

double messenger(const vector<double>& x){
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
			obj);


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
			obj);


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
			obj);


	//Free temporary memory for MGA_DSM
	for(int i = 0; i < 6; i++) {
		delete[] problem.r[i];
		delete[] problem.v[i];
	}
	problem.r.clear();
	problem.v.clear();


	return obj;
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
			obj);


	//Free temporary memory for MGA_DSM
	for(int i = 0; i < 3; i++) {
		delete[] problem.r[i];
		delete[] problem.v[i];
	}
	problem.r.clear();
	problem.v.clear();

	return obj;
}

