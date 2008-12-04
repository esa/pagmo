/*
 *  GOProblems.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/17/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include "classicobjfuns.h"
#include <math.h>


using namespace std;

double testfunction (const vector<double>& x){
	double finalvalue=0;
	for (unsigned int i=0; i<x.size(); i++){
		finalvalue += x[i];
	}
	return finalvalue;
}

double rastrigin (const vector<double>& x){
	double omega = 2.0 * M_PI;
	double value=0;
	int n = x.size();

	for (int i=0; i<n; i++){
		value += x[i]*x[i] - 10.0 * cos(omega*x[i]);
	}
	return (10.0*n + value);
}

double schwefel (const vector<double>& x){
	int n = x.size();
	double value=0;

	for (int i=0; i<n; i++){
		value += x[i] * sin(sqrt(fabs(x[i])));
		}
		return (418.9829 * n - value);
}

double ackley (const vector<double>& x){
	int n = x.size();
	double omega = 2.0 * M_PI;
	double s1=0.0, s2=0.0;
	double nepero=exp(1.0);

	for (int i=0; i<n; i++){
		s1 += x[i]*x[i];
		s2 += cos(omega*x[i]);
	}
	return -20*exp(-0.2 * sqrt(1.0/n * s1))-exp(1.0/n*s2)+ 20 + nepero;
}

double lennardjones(const vector <double>& x){
	int n = x.size();
	int atoms = (n + 6) / 3;
	vector<double> dummy(3);
	vector< vector<double> > r(atoms,dummy);				//double dimensional vector containing the atom positions
	double V = 0;											//LJ potential
	double sixth,dist;

	//we transform the decision vector x in atoms positions r
	r[0][0] = r[0][1] = r[0][2] = 0;						//x1,y1,z1 fixed
	r[1][0] = r[1][1] = 0;									//x2,y2    fixed
	r[1][2] = x[0];											//z2	   is a variable
	r[2][0] = 0;											//x3	   fixed
	r[2][1] = x[1];											//y3       is a variable
	r[2][2] = x[2];											//z3	   is a variable

	for (int i=3; i<atoms; i++){
		r[i][0] = x[3*(i-2)];
		r[i][1] = x[3*(i-2)+1];
		r[i][2] = x[3*(i-2)+2];
	}

	//We evaluate the potential
	for ( int i=0; i<(atoms-1); i++ ) {
		for ( int j=(i+1); j<atoms; j++ ) {
			dist = pow(r[i][0]-r[j][0],2) +  pow(r[i][1]-r[j][1],2) +  pow(r[i][2]-r[j][2],2);  //rij^2
			if ( dist == 0.0 ) {
				return 1e+20;	//penalty
			}
			else {
				sixth = pow(dist,-3);	//rij^-6
				V += ( pow(sixth,2) - sixth );
			}
		}
	}
	return 4 * V;
}

double rosenbrock (const vector<double>& x){
	int n = x.size();
	double value=0.0;

	for (int i=0; i<n-1; i++){
		value += 100 * (x[i]*x[i] -x[i+1])*(x[i]*x[i] -x[i+1]) + (x[i]-1)*(x[i]-1);
	}
	return value;
}

double levy(const vector<double>& x){
	int n = x.size();
	double isum = 0.0;
	double jsum = 0.0;
	double ret;
	int i, j;

	for ( j=0; j<n; j+=2 ) {
		for ( i=1; i<=5; i++ ) {
			isum += (double)(i) * cos((double)(i-1)*x[j] + (double)(i));
			jsum += (double)(i) * cos((double)(i+1)*x[j+1] + (double)(i));
		}
	}

	ret = isum*jsum;
	for ( j=0; j<n; j+=2 )
		ret += pow(x[j]+1.42513,2) + pow(x[j+1]+0.80032,2);

	return ret;
}
