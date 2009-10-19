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

// 09/2009: initial version by Francesco Biscani.

#include <string>
#include <vector>

#include "GOproblem.h"
#include "example_problem.h"

// Initialise the problem using the constructor of the base class: GOProblem(int n) will initialise a problem
// with dimension n.
exampleProb::exampleProb():GOProblem(1) {}

// Objective function of the problem. Will take a vector of values as input and will return its "fitness".
double exampleProb::objfun_(const std::vector<double> &x) const
{
	return (x[0] * x[0]);
}

// Cloning method: this will return a pointer to a new instance of the problem, copied from the current one.
// This method will look the same in virtually all problems.
exampleProb *exampleProb::clone() const
{
	return new exampleProb(*this);
}

// This method is supposed to return a description of the problem.
std::string exampleProb::id_object() const
{
	return std::string("Example problem: minimization of y = x*x.");
}