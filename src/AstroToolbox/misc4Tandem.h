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


