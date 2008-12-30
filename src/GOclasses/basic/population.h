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

#include "../../exceptions.h"
#include "GOproblem.h"
#include "constants.h"
#include "individual.h"
#include "rng.h"

class Population{
	typedef std::vector<Individual> container_type;
public:
	// The following is needed to mimic the behaviour of a std::vector, hence easing the interfacing with Python.
	typedef container_type::value_type value_type;
	typedef container_type::size_type size_type;
	typedef container_type::difference_type difference_type;
	typedef container_type::iterator iterator;
	typedef container_type::const_iterator const_iterator;
	template <class InputIterator>
	Population(InputIterator f, InputIterator l):pop(f,l) {}
	iterator begin() {return pop.begin();}
	const_iterator begin() const {return pop.begin();}
	iterator end() {return pop.end();}
	const_iterator end() const {return pop.end();}
	size_t size () const;
	iterator erase(iterator pos) {return pop.erase(pos);}
	iterator erase(iterator first, iterator last) {return pop.erase(first,last);}
	iterator insert(iterator pos, const Individual &x) {return pop.insert(pos,x);}
	template <class InputIterator>
	void insert(iterator pos, InputIterator f, InputIterator l) {return pop.insert(pos,f,l);}
	void push_back(const Individual &);
	//Methods
	Population() {}
	Population(GOProblem &, int);
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
	void createRandomPopulation(GOProblem &, int);
	template <class Functor>
	const Individual &extract_most() const {
		if (pop.empty()) {
			throw index_error("Population is empty.");
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
		const size_t size = picks.size();
		for (size_t i = 0; i < size; ++i) {
			if (Forced || deme[i].getFitness() < pop[picks[i]].getFitness()) {
				pop[picks[i]] = deme[i];
			}
		}
	}
	std::vector<Individual> pop;
};

std::ostream &operator<<(std::ostream &, const Population &);

#endif
