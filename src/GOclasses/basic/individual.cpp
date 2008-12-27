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

	Individual::Individual(const std::vector<double> &LB, const std::vector<double> &UB, rng_double &):
		x(LB.size()),v(LB.size()),fitness(0) {
		// Fill a new random chromosome and velocity vector.
		const size_t size = LB.size();
		for (size_t i = 0; i < size; ++i){
			x[i] = LB[i] + drng() * (UB[i] - LB[i]);
			v[i] = drng() * (UB[i] - LB[i]);
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

