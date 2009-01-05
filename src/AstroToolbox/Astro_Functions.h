// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //

#ifndef ASTRO_FUNCTIONS_H
#define ASTRO_FUNCTIONS_H

// Conversion from Mean Anomaly to Eccentric Anomaly via Kepler's equation
double Mean2Eccentric (const double &, const double &);

void Conversion(const double*, double*, double*, const double &);

double norm(const double*, const double*);

double norm2(const double*);

void vett(const double*, const double*, double*);

double tofabn(const double&, const double&, const double&);

void vers(const double*, double*);

double x2tof(const double&, const double&, const double&, const int &);

#endif



