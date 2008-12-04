// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //

/*
%Inputs:
%           r0:    column vector for the position (mu=1)
%           v0:    column vector for the velocity (mu=1)
%           rtarget: distance to be reached
%
%Outputs:
%           t:     time taken to reach a given distance
%
%Comments:  everything works in non dimensional units
*/

#include "time2distance.h"
#include <complex>

double time2distance(const double *r0, const double *v0, double rtarget)
{
	double temp = 0.0;
	double E[6];
	double r0norm = norm2(r0);
	double a, e, E0, p, ni, Et;
	int i;

	if (r0norm < rtarget)
	{
		for (i=0; i<3; i++)
			temp += r0[i]*v0[i];

		IC2par(r0,v0,1,E);
		a = E[0];
		e = E[1];
		E0 = E[5];
		p = a * (1-e*e);
		// If the solution is an ellipse
		if (e<1)
		{
			double ra = a * (1+e);
			if (rtarget>ra)
				return -1; // NaN;
            else // we find the anomaly where the target distance is reached
			{
				ni = acos((p/rtarget-1)/e);         //in 0-pi
				Et = 2*atan(sqrt((1-e)/(1+e))*tan(ni/2)); // algebraic kepler's equation

				if (temp>0)
					return sqrt(pow(a,3))*(Et-e*sin(Et)-E0 + e*sin(E0));
				else
				{
					E0 = -E0;
					return sqrt(pow(a,3))*(Et-e*sin(Et)+E0 - e*sin(E0));
				}
			}
		}
		else // the solution is a hyperbolae
		{
			ni = acos((p/rtarget-1)/e);         // in 0-pi
			Et = 2*atan(sqrt((e-1)/(e+1))*tan(ni/2)); // algebraic equivalent of kepler's equation in terms of the Gudermannian

			if (temp>0) // out==1
				return sqrt(pow((-a),3))*(e*tan(Et)-log(tan(Et/2+M_PI/4))-e*tan(E0)+log(tan(E0/2+M_PI/4)));
			else
			{
				E0 = -E0;
				return sqrt(pow((-a),3))*(e*tan(Et)-log(tan(Et/2+M_PI/4))+e*tan(E0)-log(tan(E0/2+M_PI/4)));
			}
		}
	}
	else
			return 12;
}
