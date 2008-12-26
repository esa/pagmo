/*
 *  PSO.h
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 Àdvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef PSO_H
#define PSO_H

#include <cmath>
#include <vector>

#include "population.h"
#include "rng.h"

class PSOalgorithm{
public:

Population evolve(Population deme, GOProblem& problem);

void initPSO(int generationsInit,
			 int SolDimInit,
			 double omegaInit,
			 double eta1Init,
			 double eta2Init,
			 double vcoeffInit,
			 uint32_t randomSeed);

private:
	int generations;
	int SolDim;
	double omega;
	double eta1;
	double eta2;
	double vcoeff;
	rng_double_type drng;
};

#endif
