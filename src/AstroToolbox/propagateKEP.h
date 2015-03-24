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

#ifndef PROPAGATEKEP_H
#define PROPAGATEKEP_H

#include "../config.h"

void __PAGMO_VISIBLE_FUNC propagateKEP(const double *, const double *, const double &, const double &,
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




