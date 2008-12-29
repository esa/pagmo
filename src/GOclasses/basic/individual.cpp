/*
 *  individual.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 Àdvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include <vector>

#include "individual.h"
#include "rng.h"

	Individual::Individual(GOProblem &problem):
		x(problem.getLB().size()),v(problem.getLB().size()),fitness(0) {
		init(problem);
		// Evaluation of fitness.
		fitness = problem.objfun(x);
	}

	// Resize the individual to the size of the problem and fill it with random values.
	void Individual::init(GOProblem &problem) {
		static_rng_double drng;
		// Store local references.
		const std::vector<double> &LB = problem.getLB(), &UB = problem.getUB();
		// Match problem's size.
		const size_t size = LB.size();
		x.resize(size);
		v.resize(size);
		// Fill a new random chromosome and velocity vector.
		for (size_t i = 0; i < size; ++i){
			x[i] = LB[i] + drng() * (UB[i] - LB[i]);
			v[i] = drng() * (UB[i] - LB[i]);
		}
	}

	std::ostream &operator<<(std::ostream &s, const Individual &ind) {
		for (size_t i = 0; i < ind.x.size(); ++i) {
			s << ind.x[i] << " ";
		}
		return s;
	}

