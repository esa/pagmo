// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //


#ifndef miscTandem_H
#define miscTandem_H

int xant(const double(&X)[5] , const double &x);
int yant(const double(&Y)[15] , const double &y);
int xantA5(const double(&X)[9] , const double &x);
int yantA5(const double(&Y)[13] , const double &y);
double interp2SF(const double (&X)[5] ,  double(&Y)[15] , const double &VINF, const double &declination);
double interp2A5(const double (&X)[5] ,  double(&Y)[15] , const double &VINF, const double &declination);
double SoyuzFregat (const double &VINF, const double &declination) ;
double Atlas501 (const double &VINF, const double &declination) ;
void ecl2equ (double (&ecl)[3],double (&equ)[3]);

#endif


