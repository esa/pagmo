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


#include "Astro_Functions.h"

void propagateKEP(const double *, const double *, double, double,
				  double *, double *);

void IC2par(const double*, const double*, double, double*);

void par2IC(const double*, double, double*, double*);

// Returns the cross product of the vectors X and Y.
// That is, z = X x Y.  X and Y must be 3 element
// vectors.
void cross(const double*, const double*, double*);

#endif




