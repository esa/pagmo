/*
 *  DE.h
 *  SeGMO
 *
 *  Created by Dario Izzo on 5/18/08.
 *  Copyright 2008 Àdvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef DE_H
#define DE_H

#include <vector>
#include <cmath>

#include "constants.h"
#include "population.h"
#include "rng.h"

class DEalgorithm{
public:

	Population evolve(Population deme, GOProblem& problem);

	void initDE(int generationsInit,
				 int SolDimInit,
				 double FInit,
				 double CRInit,
				 int strategyInit,
				 uint32_t randomSeed);

private:
	int generations;
	int SolDim;
	double F;
	double CR;
	int strategy;
	rng_double drng;
};

#endif
