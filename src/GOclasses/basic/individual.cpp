/*
 *  individual.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 Àdvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include <algorithm>
#include <iostream>
#include <vector>

#include "individual.h"
#include "rng.h"

	Individual::Individual(const GOProblem &problem):
		x(problem.getLB().size()),v(problem.getLB().size()),fitness(0) {
		init(problem);
	}

	Individual::Individual(const std::vector<double> &x_, const std::vector<double> &v_, const double &fitness_):
		x(x_),v(v_),fitness(fitness_) {
		// In case the input vectors have different sizes, make sure both x and v have the max size.
		const size_t max_size = std::max<size_t>(x.size(),v.size());
		x.resize(max_size);
		v.resize(max_size);
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

	std::ostream &operator<<(std::ostream &s, const Individual &ind) {
		for (size_t i = 0; i < ind.x.size(); ++i) {
			s << ind.x[i] << " ";
		}
		return s;
	}

