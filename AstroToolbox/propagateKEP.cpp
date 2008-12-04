// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //
//
// File: propageteKEP.cpp
//

#include "propagateKEP.h"

/*
 Origin: MATLAB code programmed by Dario Izzo (ESA/ACT)

 C++ version by Tamas Vinko (ESA/ACT)

 Inputs:
           r0:    column vector for the non dimensional position
           v0:    column vector for the non dimensional velocity
           t:     non dimensional time

 Outputs:
           r:    column vector for the non dimensional position
           v:    column vector for the non dimensional velocity

 Comments:  The function works in non dimensional units, it takes an
 initial condition and it propagates it as in a kepler motion analytically.
*/

void propagateKEP(const double *r0_in, const double *v0_in, double t, double mu,
				  double *r, double *v)
{

/*
   The matrix DD will be almost always the unit matrix, except for orbits
   with little inclination in which cases a rotation is performed so that
   par2IC is always defined
*/

	double DD[9] = {1, 0, 0,
		            0, 1, 0,
					0, 0, 1};
	double h[3];
	double ih[3] = {0,0,0};
	double temp1[3] = {0,0,0}, temp2[3] = {0,0,0};
	double E[6];
	double normh, M, M0;
	double r0[3], v0[3];

	int i;

	for (i=0; i<3; i++)
	{
			r0[i] = r0_in[i];
			v0[i] = v0_in[i];
	}

	vett(r0, v0, h);

	normh=norm2(h);

	for (i=0; i<3; i++)
		ih[i] = h[i]/normh;

	if (fabs(fabs(ih[2])-1.0) < 1e-3)     // the abs is needed in cases in which the orbit is retrograde,
	{                                 // that would held ih=[0,0,-1]!!
		DD[0] = 1; DD[1] =  0;  DD[2] = 0;
		DD[3] = 0; DD[4] =  0;  DD[5] = 1;
		DD[6] = 0; DD[7] = -1; DD[8] = 0;

		// Random rotation matrix that make the Euler angles well defined for the case
		// For orbits with little inclination another ref. frame is used.

		for (int i=0; i<3; i++)
		{
			temp1[0] += DD[i]*r0[i];
			temp1[1] += DD[i+3]*r0[i];
			temp1[2] += DD[i+6]*r0[i];
			temp2[0] += DD[i]*v0[i];
			temp2[1] += DD[i+3]*v0[i];
			temp2[2] += DD[i+6]*v0[i];
		}
		for (int i=0; i<3; i++)
		{
			r0[i] = temp1[i];
			temp1[i] = 0.0;
			v0[i] = temp2[i];
			temp2[i] = 0.0;
		}
        // for practical reason we take the transpose of the matrix DD here (will be used at the end of the function)
		DD[0] = 1; DD[1] =  0;  DD[2] = 0;
		DD[3] = 0; DD[4] =  0;  DD[5] = -1;
		DD[6] = 0; DD[7] =  1;  DD[8] = 0;
	}

	IC2par(r0, v0, mu, E);
	if (E[1] < 1.0)
	{
		M0 = E[5] - E[1]*sin(E[5]);
		M=M0+sqrt(mu/pow(E[0],3))*t;
	}
	else
	{
		M0 = E[1]*tan(E[5]) - log(tan(0.5*E[5] + 0.25*M_PI));
		M=M0+sqrt(mu/pow(-E[0],3))*t;
	}

	E[5]=Mean2Eccentric(M, E[1]);
	par2IC(E, mu, r, v);


	for (int j=0; j<3; j++)
	{
		temp1[0] += DD[j]*r[j];
		temp1[1] += DD[j+3]*r[j];
		temp1[2] += DD[j+6]*r[j];
		temp2[0] += DD[j]*v[j];
		temp2[1] += DD[j+3]*v[j];
		temp2[2] += DD[j+6]*v[j];
	}
	for (int i=0; i<3; i++)
	{
		r[i] = temp1[i];
		v[i] = temp2[i];
	}

	return;
}




/*
	Origin: MATLAB code programmed by Dario Izzo (ESA/ACT)

	C++ version by Tamas Vinko (ESA/ACT) 12/09/2006

	Inputs:
           r0:    column vector for the position
           v0:    column vector for the velocity

	Outputs:
           E:     Column Vectors containing the six keplerian parameters,
                  (a,e,i,OM,om,Eccentric Anomaly (or Gudermannian whenever e>1))

	Comments:  The parameters returned are, of course, referred to the same
	ref. frame in which r0,v0 are given. Units have to be consistent, and
	output angles are in radians
	The algorithm used is quite common and can be found as an example in Bate,
	Mueller and White. It goes singular for zero inclination
*/

void IC2par(const double *r0, const double *v0, double mu, double *E)
{
	double k[3];
	double h[3];
	double Dum_Vec[3];
	double n[3];
	double evett[3];

	double p = 0.0;
	double temp =0.0;
	double R0, ni;
	int i;

	vett(r0, v0, h);

	for (i=0; i<3; i++)
		p += h[i]*h[i];

	p/=mu;


	k[0] = 0; k[1] = 0; k[2] = 1;
	vett(k, h, n);


	for (i=0; i<3; i++)
		temp += pow(n[i], 2);

	temp = sqrt(temp);

	for (i=0; i<3; i++)
		n[i] /= temp;

	R0 = norm2(r0);

	vett(v0, h, Dum_Vec);

	for (i=0; i<3; i++)
		evett[i] = Dum_Vec[i]/mu - r0[i]/R0;

	double e = 0.0;
	for (i=0; i<3; i++)
		e += pow(evett[i], 2);

	E[0] = p/(1-e);
	E[1] = sqrt(e);
	e = E[1];

	E[2] = acos(h[2]/norm2(h));

	temp = 0.0;
	for (i=0; i<3; i++)
		temp+=n[i]*evett[i];

	E[4] = acos(temp/e);

	if (evett[2] < 0) E[4] = 2*M_PI - E[4];

	E[3] = acos(n[0]);
	if (n[1] < 0) E[3] = 2*M_PI-E[3];

    temp = 0.0;
	for (i=0; i<3; i++)
		temp+=evett[i]*r0[i];

	ni = acos(temp/e/R0);  // danger, the argument could be improper.

	temp = 0.0;
	for (i=0; i<3; i++)
		temp+=r0[i]*v0[i];

	if (temp<0.0) ni = 2*M_PI - ni;

	if (e<1.0)
		E[5] = 2.0*atan(sqrt((1-e)/(1+e))*tan(ni/2.0));  // algebraic kepler's equation
	else
		E[5] =2.0*atan(sqrt((e-1)/(e+1))*tan(ni/2.0));   // algebraic equivalent of kepler's equation in terms of the Gudermannian
}

/*
	Origin: MATLAB code programmed by Dario Izzo (ESA/ACT)

	C++ version by Tamas Vinko (ESA/ACT)

	Usage: [r0,v0] = IC2par(E,mu)

	Outputs:
           r0:    column vector for the position
           v0:    column vector for the velocity

	Inputs:
           E:     Column Vectors containing the six keplerian parameters,
                  (a (negative for hyperbolas),e,i,OM,om,Eccentric Anomaly)
           mu:    gravitational constant

	Comments:  The parameters returned are, of course, referred to the same
	ref. frame in which r0,v0 are given. a can be given either in kms or AUs,
	but has to be consistent with mu.All the angles must be given in radians.
	This function does work for hyperbolas as well.
*/

void par2IC(const double *E, double mu, double *r0, double *v0)
{
	double a = E[0];
	double e = E[1];
	double i = E[2];
	double omg = E[3];
	double omp = E[4];
	double EA = E[5];
	double b, n, xper, yper, xdotper, ydotper;
	double R[3][3];
	double cosomg, cosomp, sinomg, sinomp, cosi, sini;
	double dNdZeta;

	// Grandezze definite nel piano dell'orbita

	if (e<1.0)
	{
		b = a*sqrt(1-e*e);
		n = sqrt(mu/(a*a*a));
		xper=a*(cos(EA)-e);
		yper=b*sin(EA);

		xdotper = -(a*n*sin(EA))/(1-e*cos(EA));
		ydotper=(b*n*cos(EA))/(1-e*cos(EA));
	}
	else
	{
		b = -a*sqrt(e*e-1);
		n = sqrt(-mu/(a*a*a));

		dNdZeta = e * (1+tan(EA)*tan(EA))-(0.5+0.5*pow(tan(0.5*EA + 0.25*M_PI),2))/tan(0.5*EA+ 0.25*M_PI);

		xper = a/cos(EA) - a*e;
		yper = b*tan(EA);

		xdotper = a*tan(EA)/cos(EA)*n/dNdZeta;
		ydotper = b/pow(cos(EA), 2)*n/dNdZeta;
	}

    // Matrice di trasformazione da perifocale a ECI

	cosomg = cos(omg);
	cosomp = cos(omp);
	sinomg = sin(omg);
	sinomp = sin(omp);
	cosi = cos(i);
	sini = sin(i);


	R[0][0]=cosomg*cosomp-sinomg*sinomp*cosi;
	R[0][1]=-cosomg*sinomp-sinomg*cosomp*cosi;
	R[0][2]=sinomg*sini;
	R[1][0]=sinomg*cosomp+cosomg*sinomp*cosi;
	R[1][1]=-sinomg*sinomp+cosomg*cosomp*cosi;
	R[1][2]=-cosomg*sini;
	R[2][0]=sinomp*sini;
	R[2][1]=cosomp*sini;
	R[2][2]=cosi;

	// Posizione nel sistema inerziale


	double temp[3] = {xper, yper, 0.0};
	double temp2[3] = {xdotper, ydotper, 0};

	for (int j = 0; j<3; j++)
	{
		r0[j] = 0.0; v0[j] = 0.0;
		for (int k = 0; k<3; k++)
		{
			r0[j]+=R[j][k]*temp[k];
			v0[j]+=R[j][k]*temp2[k];
		}
	}
	return;
}

// Returns the cross product of the vectors X and Y.
// That is, z = X x Y.  X and Y must be 3 element
// vectors.
void cross(const double *x, const double *y, double *z)
{
   z[0] = x[1]*y[2] - x[2]*y[1];
   z[1] = x[2]*y[0] - x[0]*y[2];
   z[2] = x[0]*y[1] - x[1]*y[0];
}
