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

#ifndef ASTRO_FUNCTIONS_H
#define ASTRO_FUNCTIONS_H

#include "../config.h"

// Conversion from Mean Anomaly to Eccentric Anomaly via Kepler's equation
double __PAGMO_VISIBLE_FUNC Mean2Eccentric (const double &, const double &);

void Conversion(const double*, double*, double*, const double &);

double norm(const double*, const double*);

double norm2(const double*);

void vett(const double*, const double*, double*);

double tofabn(const double&, const double&, const double&);

void vers(const double*, double*);

double x2tof(const double&, const double&, const double&, const int &);

#endif



