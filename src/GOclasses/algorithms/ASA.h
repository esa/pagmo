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

	Population evolve(const Individual &, GOProblem &);

	//This method initialise all the SA-AN algorithm parameters
	void initASA(int niterTotInit,
				 int niterTempInit,
				 int niterRangeInit,
				 int SolDimInit,
				 double T0Init,
				 double TcoeffInit,
				 double StartStepInit,
				 uint32_t randomSeed);

	//This method initialise the SA-AN algorithm starting and final temperature setting deafult values for
	//the StartStep, the niterTemp and the niterRange. Tcoeff is evaluated accordingly
	void initASA(int niterTotInit,
				 int SolDimInit,
				 double Ts,
				 double Tf,
				 uint32_t randomSeed);

private:
	int niterTot;
	int niterTemp;
	int niterRange;
	int SolDim;
	double T0;
	double Tcoeff;
	double StartStep;
    int niterOuter;
    rng_double drng;
};

#endif
