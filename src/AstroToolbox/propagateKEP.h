// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2006 European Space Agency                            //
// ------------------------------------------------------------------------ //

#ifndef PROPAGATEKEP_H
#define PROPAGATEKEP_H

void propagateKEP(const double *, const double *, const double &, const double &,
				  double *, double *);

void IC2par(const double*, const double*, const double &, double*);

void par2IC(const double*, const double &, double*, double*);

// Returns the cross product of the vectors X and Y.
// That is, z = X x Y.  X and Y must be 3 element
// vectors.
inline void cross(const double *x, const double *y, double *z)
{
   z[0] = x[1]*y[2] - x[2]*y[1];
   z[1] = x[2]*y[0] - x[0]*y[2];
   z[2] = x[0]*y[1] - x[1]*y[0];
}

#endif




