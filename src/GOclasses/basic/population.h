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

#include "GOproblem.h"
#include "constants.h"
#include "individual.h"
#include "rng.h"

class Population{
public:
    //Methods
	void createRandomPopulation(const std::vector<double> &, const std::vector<double> &, int, rng_double &);
	void evaluatePopulation(GOProblem &);
	void addIndividual(const Individual &);
	void substituteIndividual(const Individual &, int);
	double evaluateMean() const;
	double evaluateStd() const;
	unsigned int size () const;
	Individual extractBestIndividual() const;
	Individual extractWorstIndividual() const;
	Population extractRandomDeme(int, std::vector<int> &, rng_double &);
	void insertDeme(const Population &, const std::vector<int> &);
	void insertBestInDeme(const Population &, const std::vector<int> &);
	void insertDemeForced(const Population &, const std::vector<int> &);

	//Operators
	Individual &operator[](int);
	const Individual &operator[](int) const;
	void operator=(const Population &newpop);
	void operator=(const Individual &x);

	//logging function
	friend std::ostream& operator<<(std::ostream& s, Population& pop);
private:
	std::vector<Individual> pop;
};//class Population

std::ostream& operator<<(std::ostream& s, Population& pop);
#endif
