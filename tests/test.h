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

#ifndef PAGMO_TEST_H
#define PAGMO_TEST_H

#include "../src/pagmo.h"

#define PRINT_VEC(x) do{ std::cout<<#x "[] = "; for(unsigned int __i=0;__i<(x).size();__i++) std::cout<<(x)[__i]<<" "; std::cout<<std::endl; } while(0);

bool is_eq(double d1, double d2, double eps=10e-9)
{
	return (fabs(d1 - d2) <= eps);
}

bool is_eq_vector(const pagmo::fitness_vector &f1, const pagmo::fitness_vector &f2, double eps=10e-9)
{
	if(f1.size() != f2.size()) return false;
	for(pagmo::fitness_vector::size_type i = 0; i < f1.size(); i++) {
		if (!is_eq(f1[i], f2[i], eps)) return false;
	}
	return true;
}

#endif
