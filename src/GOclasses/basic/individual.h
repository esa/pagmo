/*
 *  individual.h
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 Advanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <iostream>
#include <vector>

#include "GOproblem.h"
#include "rng.h"

class Individual{

public:
	Individual(const std::vector<double> &, const std::vector<double> &, const double &);
	Individual(const GOProblem &);
	double getFitness() const {return fitness;}
	const std::vector<double> &getDecisionVector() const {return x;}
	const std::vector<double> &getVelocity() const {return v;}
	friend std::ostream &operator<<(std::ostream &, const Individual &);
private:
	void init(const GOProblem &);
	std::vector<double> x;  //this is the "chromosome" or "decision vector"
	std::vector<double> v;  //this is the "velocity" or "heading" of each individual
	double fitness;
};

std::ostream &operator<<(std::ostream &, const Individual &);

#endif
