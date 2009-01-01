/*
 *  ASA.h
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef ASA_H
#define ASA_H

#include <vector>
#include <cmath>

#include "GOproblem.h"
#include "constants.h"
#include "individual.h"
#include "population.h"
#include "rng.h"


class ASAalgorithm{
public:
	//This method initialise the SA-AN algorithm starting and final temperature setting deafult values for
	//the StartStep, the niterTemp and the niterRange. Tcoeff is evaluated accordingly.
	ASAalgorithm(int, const GOProblem &, const double &, const double &);
	Population evolve(const Population &, GOProblem &);
private:
	size_t niterTot;
	size_t niterTemp;
	size_t niterRange;
	const size_t SolDim;
	double T0;
	double Tcoeff;
	double StartStep;
	size_t niterOuter;
	rng_double drng;
};

#endif
