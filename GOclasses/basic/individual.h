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
#include "constants.h"
#include "rng.h"

class Individual{

public:
    //methods
	Individual(const std::vector<double> &x_, const std::vector<double> &v_, const double &fitness_):x(x_),v(v_),fitness(fitness_) {}
	Individual(const std::vector<double> &, const std::vector<double> &, rng_double_type &);
	double evaluateFitness(GOProblem &);
	double getFitness() const;
	const std::vector<double> &getDecisionVector() const;
	void setDecisionVector(const std::vector<double> &);
	const std::vector<double> &getVelocity() const;
	void setVelocity(const std::vector<double> &);
	void resetVelocity(const std::vector<double> &, const std::vector<double> &, rng_double_type &);

	//operators
	double &operator[](int);
	const double &operator[](int) const;
	void operator=(const Individual &);

	//logging function
	friend std::ostream& operator<<(std::ostream &, const Individual &);
private:
	std::vector<double> x;	//this is the "chromosome" or "decision vector"
	std::vector<double> v;  //this is the "velocity" or "heading" of each individual
	double fitness;
};// class Individual

std::ostream& operator<<(std::ostream &, const Individual &);
#endif
