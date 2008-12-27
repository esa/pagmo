/*
 *  individual.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 Àdvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include "individual.h"
#include "rng.h"

	mt_rng_double Individual::drng;

	Individual::Individual(GOProblem &problem):
		x(problem.getLB().size()),v(problem.getLB().size()),fitness(0) {
		// Fill a new random chromosome and velocity vector.
		const size_t size = problem.getLB().size();
		for (size_t i = 0; i < size; ++i){
			x[i] = problem.getLB()[i] + drng() * (problem.getUB()[i] - problem.getLB()[i]);
			v[i] = drng() * (problem.getUB()[i] - problem.getLB()[i]);
		}
	}

	const double &Individual::evaluateFitness(GOProblem &problem) {
		fitness = problem.objfun(x);
		return fitness;
	};

	std::ostream &operator<<(std::ostream &s, const Individual &ind) {
		for (size_t i = 0; i < ind.x.size(); ++i) {
			s << ind.x[i] << " ";
		}
		return s;
	}

