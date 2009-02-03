/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

// 04/06/08 Created by Dario Izzo.

#ifndef PAGMO_GOPROBLEM_H
#define PAGMO_GOPROBLEM_H

#include <string>
#include <typeinfo>
#include <vector>

#include "../../config.h"

class __PAGMO_VISIBLE GOProblem {
public:
	// Virtual destructor - required because the class contains a pure virtual member function
	virtual ~GOProblem() {};
	// Bounds getters and setters via reference
	const std::vector<double> &getLB() const {return LB;}
	const std::vector<double> &getUB() const {return UB;}
	// Dimension getter
	size_t getDimension() const {return LB.size();}
	// The objective function - must be implemented in subclasses
	virtual double objfun(const std::vector<double> &) const = 0;
	virtual GOProblem *clone() const = 0;
	std::string id_name() const {return typeid(*this).name();}
protected:
	// Constructor with array bounds initialisers
	GOProblem(const size_t &d, const double *l, const double *u):LB(l, l + d),UB(u, u + d) {};
	// Constructor with vectors initialisers
	GOProblem(const std::vector<double>& l, const std::vector<double>& u):LB(l),UB(u) {};
	//Default Constructor. Necessary as some derived classes need to call it (i.e. LJ)
	GOProblem() {};
	//These need to be protected and cannot be private and const as some problems need to deifne the LB and UB at run time (i.e. LJ)
	std::vector<double> LB;
	std::vector<double> UB;
};

#endif
