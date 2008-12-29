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
	Population() {}
	Population(GOProblem &, int);
	void push_back(const Individual &);
	double evaluateMean() const;
	double evaluateStd() const;
	size_t size () const;
	const Individual &extractBestIndividual() const;
	const Individual &extractWorstIndividual() const;
	Population extractRandomDeme(int, std::vector<int> &, rng_double &);
	void insertDeme(const Population &, const std::vector<int> &);
	void insertBestInDeme(const Population &, const std::vector<int> &);
	void insertDemeForced(const Population &, const std::vector<int> &);

	//Operators
	Individual &operator[](const size_t &);
	const Individual &operator[](const size_t &) const;

private:
	void createRandomPopulation(GOProblem &, int);
	template <class Functor>
	const Individual &extract_most() const {
		const Functor func = Functor();
		double f = pop[0].getFitness();
		size_t index = 0;
		const size_t size = pop.size();
		for (size_t i = 1; i < size; ++i) {
			if(func(pop[i].getFitness(),f)) {
				index = i;
				f = pop[i].getFitness();
			}
		}
		return pop[index];
	}
	std::vector<Individual> pop;
};//class Population

std::ostream &operator<<(std::ostream &, const Population &);
#endif
