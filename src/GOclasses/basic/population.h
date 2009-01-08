/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

// 16/05/08 Created by Dario Izzo.

#ifndef POPULATION_H
#define POPULATION_H

#include <boost/scoped_ptr.hpp>
#include <vector>

#include "../../config.h"
#include "../../exceptions.h"
#include "GOproblem.h"
#include "constants.h"
#include "individual.h"
#include "py_container_utils.h"
#include "rng.h"

class __PAGMO_VISIBLE Population: public py_container_utils<Population> {
public:
	//Methods
	Population(const GOProblem &);
	Population(const GOProblem &, int);
	Population(const Population &);
	Population &operator=(const Population &);
	Individual &operator[](int);
	const Individual &operator[](int) const;
	void push_back(const Individual &);
	void insert(int, const Individual &);
	void erase(int);
	size_t size() const;
	const GOProblem &problem() const {return *m_problem;}
	double evaluateMean() const;
	double evaluateStd() const;
	Individual extractBestIndividual() const;
	Individual extractWorstIndividual() const;
	// TODO: review this API.
	Population extractRandomDeme(int, std::vector<size_t> &);
	void insertDeme(const Population &, const std::vector<size_t> &);
	void insertBestInDeme(const Population &, const std::vector<size_t> &);
	void insertDemeForced(const Population &, const std::vector<size_t> &);
private:
	void check_individual(const Individual &) const;
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
			pagmo_throw(index_error,"mismatch between deme size and picks size while inserting deme");
		}
		for (size_t i = 0; i < picks_size; ++i) {
			if (picks[i] >= pop_size) {
				pagmo_throw(index_error,"pick value exceeds population's size while inserting deme");
			}
			if (Forced || deme[i].getFitness() < pop[picks[i]].getFitness()) {
				pop[picks[i]] = deme[i];
			}
		}
	}
	std::vector<Individual>				pop;
	boost::scoped_ptr<const GOProblem>	m_problem;
};

inline std::ostream &operator<<(std::ostream &s, const Population &pop) {
	s << "Problem type: '" << pop.problem().id_name() << "'\n";
	for (size_t i = 0; i < pop.size(); ++i) {
		s << "Individual #" << i << ": " << pop[i].getFitness() << " " << pop[i] << std::endl;
	}
	return s;
}

#endif
