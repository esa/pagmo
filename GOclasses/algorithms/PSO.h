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

#include "population.h"
#include <vector>
#include <math.h>
#include "PkRandom.h"

class PSOalgorithm{
public:

Population evolve(Population deme, GOProblem& problem);

void initPSO(int generationsInit,
			 int SolDimInit,
			 double omegaInit,
			 double eta1Init,
			 double eta2Init,
			 double vcoeffInit,
			 unsigned long randomSeed);

private:
	int generations;
	int SolDim;
	double omega;
	double eta1;
	double eta2;
	double vcoeff;
	Pk::Random32 rng;
};

#endif
