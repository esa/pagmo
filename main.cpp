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
 *   the Free Software Foundation; either version 3 of the License, or       *
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

#include <iostream>
#include <cmath>
#include <gal_odeint.h>

using namespace std;


void dy (double t,double y[],double dy[], int* param);

int main()
{
        double y[6]; 			//Contains the spacecraft state
	int nvar = 6;	  		//This is the number of equations
	double t0 = 0;	  		//Starting integration time
	double tf = 2 * M_PI;		//End integration time
	double eps = 1e-14;		//Accuracy
	double h1 = 1e-2;		//First guess for the step
	double hmin = 1e-19;		//Minimum allowed step size
	
	int* param;			//Parameters
	int retval;			//RetVal
	
	double thrust[3];
	thrust[0] = 1e-4;
	thrust[1] = 0;
	thrust[2] = 0;
	
	param = (int*) thrust;
	
	//Initialise initial conditions
	y[0] = 0.1;
	y[1] = 0;
	y[2] = 0;
	y[3] = 0;
	y[4] = sqrt(1/y[0]);
	y[5] = 0;
	
	//Set integration tuime to one period
	tf *= pow(y[0],3.0/2.0);
	
	cout << "Test for the GAL library ODE integrators: " << endl;
	cout << "Initial conditions are: " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << endl;
	
	//Integrate
	retval = gal_rkf(y,nvar,t0,tf,eps,h1,hmin,dy,gal_rkfs78,param);
	
	cout << "Return value: " << retval << endl;
	cout << "Final conditions are: " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << endl;
        return 0;
}

void dy (double t,double y[],double dy[], int* param){

  double* thrust;
  thrust = (double*) param;
  
  double r = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
  double r3 = r*r*r;
  dy[0] = y[3];
  dy[1] = y[4];
  dy[2] = y[5];
  dy[3] = - y[0] / r3 + thrust[0];
  dy[4] = - y[1] / r3 + thrust[1];
  dy[5] = - y[2] / r3 + thrust[2];
}