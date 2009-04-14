/*
 *  SGA.h
 *  Simple Genetic Algorithm
 *
 *  Created by Dario Izzo on 10/5/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SGA_H
#define SGA_H

#include <cmath>
#include <vector>

#include "GOproblem.h"
#include "population.h"
#include "rng.h"

class SGAalgorithm{
public:

Population evolve(Population deme, GOProblem& problem);

void initSGA(int generationsInit,
			 int SolDimInit,
			 double CRInit,
			 double MInit,
			 int insert_bestInit,
			 uint32_t randomSeed
			 );
	
	//virtual std::string id_object() const {return id_name(); }

private:
	int generations;
	int SolDim;
	double CR;		//crossover
	double M;		//mutation
	int insert_best;
	rng_double drng;
	rng_uint32 rng;
};

#endif
