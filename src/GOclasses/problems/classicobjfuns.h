/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

// 17/05/08 Created by Dario Izzo.

#ifndef PAGMO_CLASSICOBJFUNS_H
#define PAGMO_CLASSICOBJFUNS_H

#include <vector>

#include "../../config.h"

//NOTE: the functions here have passing by reference + const as they are called a lot of time during execution and thus
//it is worth trying to save time by avoiding to make a copy of the variable passed

double __PAGMO_VISIBLE_FUNC testfunction (const std::vector<double>& x);
double __PAGMO_VISIBLE_FUNC rastriginf (const std::vector<double>& x);
double __PAGMO_VISIBLE_FUNC schwefelf (const std::vector<double>& x);
double __PAGMO_VISIBLE_FUNC ackleyf (const std::vector<double>& x);
double __PAGMO_VISIBLE_FUNC rosenbrockf (const std::vector<double>& x);
double __PAGMO_VISIBLE_FUNC lennardjonesf (const std::vector<double>& x);
double __PAGMO_VISIBLE_FUNC levyf (const std::vector<double>& x);
double __PAGMO_VISIBLE_FUNC griewankf (const std::vector<double>& x);

#endif
