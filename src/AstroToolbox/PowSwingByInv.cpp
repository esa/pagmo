/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#include "PowSwingByInv.h"
#include <math.h>

void PowSwingByInv (const double Vin,const double Vout,const double alpha,
					double &DV,double &rp)
{
	const int maxiter = 30;
	int i = 0;
	double err = 1.0;
	double f,df;                    // function and its derivative
	double rp_new;
	const double tolerance = 1e-8;

	double aIN  = 1.0/pow(Vin,2);     // semimajor axis of the incoming hyperbola
	double aOUT = 1.0/pow(Vout,2);    // semimajor axis of the incoming hyperbola

	rp = 1.0;
	while ((err > tolerance)&&(i < maxiter))
	{
		i++;
		f = asin(aIN/(aIN + rp)) + asin(aOUT/(aOUT + rp)) - alpha;
		df = -aIN/sqrt((rp + 2 * aIN) * rp)/(aIN+rp) -
			 aOUT/sqrt((rp + 2 * aOUT) * rp)/(aOUT+rp);
		rp_new = rp - f/df;
		if (rp_new > 0)
		{
			err = fabs(rp_new - rp);
			rp = rp_new;
		}
		else
			rp /= 2.0;
	}

	// Evaluation of the DV
	DV = fabs (sqrt(Vout*Vout + (2.0/rp)) - sqrt(Vin*Vin + (2.0/rp)));
}
