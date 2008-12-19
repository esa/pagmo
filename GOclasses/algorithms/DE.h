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

#include "population.h"
#include <vector>
#include <math.h>
#include "constants.h"
#include "PkRandom.h"

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
	rng_double_type drng;
};

#endif
