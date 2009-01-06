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

// 16/05/08 Created by Dario Izzo.

#include <algorithm>
#include <iostream>
#include <vector>

#include "../../exceptions.h"
#include "individual.h"
#include "rng.h"

	Individual::Individual(const GOProblem &problem):
		x(problem.getLB().size()),v(problem.getLB().size()),fitness(0) {
		init(problem);
	}

	Individual::Individual(const std::vector<double> &x_, const std::vector<double> &v_, const double &fitness_):
		x(x_),v(v_),fitness(fitness_) {
		if (x.size() != v.size()) {
			pagmo_throw(value_error,"while constructing individual, size mismatch between decision vector and velocity vector");
		}
	}

	// Resize the individual to the size of the problem and fill it with random values.
	void Individual::init(const GOProblem &problem) {
		static_rng_double drng;
		// Store local references.
		const std::vector<double> &LB = problem.getLB(), &UB = problem.getUB();
		const size_t size = LB.size();
		// Fill a new random chromosome and velocity vector.
		for (size_t i = 0; i < size; ++i){
			x[i] = LB[i] + drng() * (UB[i] - LB[i]);
			v[i] = drng() * (UB[i] - LB[i]);
		}
		// Evaluation of fitness.
		fitness = problem.objfun(x);
	}
