/*
 *  population.h
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef POPULATION_H
#define POPULATION_H

#include <vector>
#include "constants.h"
#include "individual.h"
#include "PkRandom.h"
#include "GOproblem.h"

class Population{
public:
    //Methods
	void createRandomPopulation(std::vector<double> LB, std::vector<double> UB, int N, Pk::Random32& rng);
	void resetVelocities(std::vector<double> LB, std::vector<double> UB, Pk::Random32& rng);
	void evaluatePopulation(GOProblem&);
	void addIndividual(Individual x);
	void substituteIndividual(const Individual x, const int n);
	double evaluateMean();
	double evaluateStd();
	unsigned int size ();
	Individual extractBestIndividual();
	Individual extractWorstIndividual();
	Population extractRandomDeme(int N, std::vector<int> &picks, Pk::Random32& rng);
	void insertDeme(Population deme, std::vector<int> picks);
	void insertBestInDeme(Population deme, std::vector<int> picks);
	void insertDemeForced(Population deme, std::vector<int> picks);

	//Operators
	Individual& operator[](int index);
	void operator=(Population newpop);
	void operator=(Individual x);

	//logging function
	friend std::ostream& operator<<(std::ostream& s, Population& pop);
private:
	std::vector<Individual> pop;
};//class Population

std::ostream& operator<<(std::ostream& s, Population& pop);
#endif
