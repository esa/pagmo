// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //

#include <cmath>
#include <valarray>
#include <vector>
#include "ZeroFinder.h"

namespace ZeroFinder {
	void Function1D::SetParameters(double a, double b)
	{
	  p1 = a;
	  p2 = b;
	}

	void Function1D_7param::SetParameters(double a, double b, double c, double d, double e, double f, double g)
	{
	  p1 = a;
	  p2 = b;
	  p3 = c;
	  p4 = d;
	  p5 = e;
	  p6 = f;
	  p7 = g;
	}
}

