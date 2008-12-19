/*
 *  MPSO.h
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef MPSO_H
#define MPSO_H

#include "population.h"
#include <vector>
#include <math.h>
#include "PkRandom.h"
#include "GOproblem.h"

class MPSOalgorithm{
public:

Population evolve(Population deme, GOProblem& problem);

void initMPSO(int generationsInit,
			 int SolDimInit,
			 double omegaInit,
			 double eta1Init,
			 double eta2Init,
			 double vcoeffInit,
			 int nswarmsInit,
			 uint32_t randomSeed);

private:
	int generations;
	int SolDim;
	double omega;
	double eta1;
	double eta2;
	double vcoeff;
	int nswarms;
	rng_double_type drng;
};

#endif
