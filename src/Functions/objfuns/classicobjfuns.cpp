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

// 17/05/2008: Initial version by Dario Izzo.

#include <cmath>

#include "../../constants.h"
#include "classicobjfuns.h"

using namespace std;

double testfunction (const vector<double>& x){
	double finalvalue=0;
	for (unsigned int i=0; i<x.size(); i++){
		finalvalue += x[i];
	}
	return finalvalue;
}

double rastriginf (const vector<double>& x){
	double omega = 2.0 * M_PI;
	double value=0;
	int n = x.size();

	for (int i=0; i<n; i++){
		value += x[i]*x[i] - 10.0 * cos(omega*x[i]);
	}
	return (10.0*n + value);
}

double schwefelf (const vector<double>& x){
	int n = x.size();
	double value=0;

	for (int i=0; i<n; i++){
		value += x[i] * sin(sqrt(fabs(x[i])));
		}
		return (418.9829 * n - value);
}

double ackleyf (const vector<double>& x){
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

//lennardjones helper function that transforms the decision vector x in atoms positions r
static double r(const int& atom, const int& coord, const vector <double>& x) {
	if(atom == 0) { //x1,y1,z1 fixed
		return 0.0;
	} else if(atom == 1) {
		if(coord < 2) { //x2,y2    fixed
			return 0.0;
		} else { //z2 is a variable
			return x[0];
		}
	} else if(atom == 2) {
		if(coord == 0) { //x3	   fixed
			return 0.0;
		} else { //y3 and x3 are variables
			return x[coord];
		}
	} else {
		return x[3 * (atom - 2) + coord];
	}
}

double lennardjonesf(const vector <double>& x){
	int n = x.size();
	int atoms = (n + 6) / 3;

        double V = 0;         //LJ potential
	double sixth, dist;
	
	//We evaluate the potential
	for ( int i=0; i<(atoms-1); i++ ) {
		for ( int j=(i+1); j<atoms; j++ ) {
			dist = pow(r(i, 0, x) - r(j, 0, x), 2) + pow(r(i, 1, x) - r(j, 1, x), 2) + pow(r(i, 2, x) - r(j, 2, x), 2);  //rij^2
			if ( dist == 0.0 ) {
				return 1e+20;	//penalty
			}
			else {
				sixth = pow(dist, -3);	//rij^-6
				V += (pow(sixth, 2) - sixth);
			}
		}
	}
	return 4 * V;
}

double rosenbrockf (const vector<double>& x){
	int n = x.size();
	double value=0.0;

	for (int i=0; i<n-1; i++){
		value += 100 * (x[i]*x[i] -x[i+1])*(x[i]*x[i] -x[i+1]) + (x[i]-1)*(x[i]-1);
	}
	return value;
}

double levyf(const vector<double>& x){
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

double griewankf (const vector<double>& x){
        int n = x.size();
        double fr=4000.0;
        double retval = 0.0;
        double p = 1.0;

        for (int i=0; i<n; i++){ retval += x[i]*x[i];}
        for (int i=0; i<n; i++){ p *= cos(x[i]/sqrt(x[i]));}
        return (retval/fr - p + 1);
}
