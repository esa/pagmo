// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //

#ifndef PL_EPH_AN_H
#define PL_EPH_AN_H

#include <vector>

using namespace std;

void Planet_Ephemerides_Analytical (const double, const int, double*, double*);

void Custom_Eph(const double, const double, const double[], double*, double*);

#endif
