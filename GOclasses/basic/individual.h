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
	void createRandomIndividual(std::vector<double> LB, std::vector<double> UB, rng_double_type &);
	double evaluateFitness(GOProblem&);
	double getFitness() const;
	void setFitness(double fitnessnew);
	std::vector<double> getDecisionVector() const;
	void setDecisionVector(std::vector<double> xnew);
	std::vector<double> getVelocity() const;
	void setVelocity(std::vector<double> xnew);
	void resetVelocity(std::vector<double> LB, std::vector<double> UB, rng_double_type &);

	//operators
	double& operator[](int index);
	void operator=(const Individual& newindividual);

	//logging function
	friend std::ostream& operator<<(std::ostream& s, Individual& x);
private:
	std::vector<double> x;	//this is the "chromosome" or "decision vector"
	std::vector<double> v;  //this is the "velocity" or "heading" of each individual
	double fitness;
};// class Individual

std::ostream& operator<<(std::ostream& s, Individual& x);
#endif
