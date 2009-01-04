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
#include "go_algorithm.h"
#include "individual.h"
#include "population.h"

class ASAalgorithm: public go_algorithm {
public:
	//This method initialise the SA-AN algorithm starting and final temperature setting deafult values for
	//the StartStep, the niterTemp and the niterRange. Tcoeff is evaluated accordingly.
	ASAalgorithm(int, const double &, const double &);
	virtual Population evolve(const Population &, const GOProblem &);
	virtual ASAalgorithm *clone() const {return new ASAalgorithm(*this);}
private:
	const size_t niterTot;
	const size_t niterTemp;
	const size_t niterRange;
	const double Ts;
	const double Tf;
	const double StartStep;
	size_t niterOuter;
};

#endif
