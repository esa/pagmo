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

#include <boost/scoped_ptr.hpp>
#include <vector>

#include "../../exceptions.h"
#include "GOproblem.h"
#include "constants.h"
#include "individual.h"
#include "rng.h"

class Population{
public:
	//Methods
	Population(const GOProblem &);
	Population(const GOProblem &, int);
	Population(const Population &);
	Population &operator=(const Population &);
	void push_back(const Individual &);
	size_t size () const;
	const GOProblem &problem() const {return *m_problem;}
	double evaluateMean() const;
	double evaluateStd() const;
	const Individual &extractBestIndividual() const;
	const Individual &extractWorstIndividual() const;
	Population extractRandomDeme(int, std::vector<size_t> &);
	void insertDeme(const Population &, const std::vector<size_t> &);
	void insertBestInDeme(const Population &, const std::vector<size_t> &);
	void insertDemeForced(const Population &, const std::vector<size_t> &);
	//Operators
	Individual &operator[](const size_t &);
	const Individual &operator[](const size_t &) const;
private:
	void createRandomPopulation(int);
	template <class Functor>
	const Individual &extract_most() const {
		if (pop.empty()) {
			pagmo_throw(index_error,"population is empty");
		}
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
	template <bool Forced>
	void ll_insert_deme(const Population &deme, const std::vector<size_t> &picks) {
		const size_t picks_size = picks.size(), pop_size = size();
		if (picks_size != deme.size()) {
			pagmo_throw (index_error,"mismatch between deme size and picks size while inserting deme");
		}
		for (size_t i = 0; i < picks_size; ++i) {
			if (picks[i] >= pop_size) {
				pagmo_throw (index_error,"pick value exceeds population's size while inserting deme");
			}
			if (Forced || deme[i].getFitness() < pop[picks[i]].getFitness()) {
				pop[picks[i]] = deme[i];
			}
		}
	}
	std::vector<Individual>			pop;
	boost::scoped_ptr<const GOProblem>	m_problem;
};

std::ostream &operator<<(std::ostream &, const Population &);

#endif
