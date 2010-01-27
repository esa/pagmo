/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <typeinfo>

#include "../../exceptions.h"
#include "../../rng.h"
#include "../problems/base.h"
#include "individual.h"
#include "population.h"

namespace pagmo
{

population::population(const problem::base &p)
		:m_problem(p.clone())
{
}

population::population(const problem::base &p, int N)
		:m_problem(p.clone())
{
	createRandompopulation(N);
}

population::population(const population &p)
		:pop(p.pop),
		m_problem(p.m_problem->clone())
{
	//copy costructor if the std::vector class creates a deep copy of the individuals
}

population &population::operator=(const population &p)
{
	if (this != &p) { //handle p = p
		if ((*m_problem) != (*p.m_problem)) {
			pagmo_throw(type_error, "problem types are not compatible while assigning population");
		}
		pop = p.pop; //deep copy
		m_problem.reset(p.m_problem->clone());
	}
	return *this;
}



const individual &population::operator[](int index) const
{
	return pop[get_ra_index(index)];
}

individual &population::operator[](int index)
{
	return pop[get_ra_index(index)];
}

void population::setIndividual(int idx, const individual &ind)
{
	pop[get_ra_index(idx)] = checked_individual(ind);
}

void population::push_back(const individual &i)
{
	pop.push_back(checked_individual(i));
}

void population::insert(int n, const individual &i)
{
	individual tmp = checked_individual(i);
	try {
		pop.insert(pop.begin() + get_ra_index(n),tmp);
	} catch (const index_error &) {
		if (n >= 0) {
			pop.insert(pop.begin(),tmp);
		} else {
			pop.insert(pop.end(),tmp);
		}
	}
}

void population::erase(int n)
{
	pop.erase(pop.begin() + get_ra_index(n));
}

size_t population::size() const
{
	return pop.size();
}



const problem::base &population::problem() const
{
	return *m_problem;
}



double population::evaluateMean() const
{
	if (pop.empty()) {
		pagmo_throw(index_error,"population is empty");
	}
	double mean = 0;
	const size_t size = pop.size();
	for (size_t i = 0; i < size; ++i) {
		mean += pop[i].get_fitness();
	}
	mean /= (double)size;
	return mean;
}

double population::evaluateStd() const
{
	if (pop.empty()) {
		pagmo_throw(index_error,"population is empty");
	}
	double Std = 0, mean = evaluateMean();
	const size_t size = pop.size();
	for (size_t i = 0; i < size; ++i) {
		Std += pow((pop[i].get_fitness()-mean),2.0);
	}
	Std = sqrt(Std/size);
	return Std;
}



const individual &population::extractBestIndividual() const
{
	return pop[extract_most_index<std::less<double> >()];
}

const individual &population::extractWorstIndividual() const
{
	return pop[extract_most_index<std::greater<double> >()];
}


void population::replace_best(const individual &ind)
{
	setIndividual(extract_most_index<std::less<double> >(),ind);
}

void population::replace_worst(const individual &ind)
{
	setIndividual(extract_most_index<std::greater<double> >(),ind);
}


struct fitness_sorter {
	bool operator()(const individual &i1, const individual &i2) const {
		return (i1.get_fitness() < i2.get_fitness());
	}
};

void population::sort()
{
	std::sort(pop.begin(),pop.end(),fitness_sorter());
}

population population::extractRandomDeme(int N, std::vector<size_t> &picks)
{
	if (N > (int)size()) {
		pagmo_throw(index_error,"cannot extract deme whose size is greater than the original population");
	}
	// Empty picks first.
	picks.clear();
	static_rng_double drng;
	population deme(*m_problem,0);
	const size_t pop_size = size();
	std::vector<size_t> PossiblePicks;
	PossiblePicks.reserve(pop_size);
	for (size_t i = 0; i < pop_size; ++i) {
		PossiblePicks.push_back(i);
	}
	for (int i = 0; i < N; ++i) {
		//we pick a random position between 0 and popsize-1
		const size_t Pick = (size_t)(drng() * PossiblePicks.size());
		//and store it
		picks.push_back(PossiblePicks[Pick]);
		//we insert the corresponding individual in the deme
		deme.push_back(pop[PossiblePicks[Pick]]);
		//and erase it from the possible picks
		PossiblePicks.erase(PossiblePicks.begin() + Pick);
	}
	return deme;
}

void population::insertDeme(const population &deme, const std::vector<size_t> &picks)
{
	ll_insert_deme<false>(deme,picks);
}

void population::insertDemeForced(const population &deme, const std::vector<size_t> &picks)
{
	ll_insert_deme<true>(deme,picks);
}

void population::insertBestInDeme(const population &deme, const std::vector<size_t> &picks)
{
	const size_t Ndeme = deme.size(), pop_size = size();
	if (picks.size() != Ndeme) {
		pagmo_throw(index_error,"mismatch between deme size and picks size while inserting best in deme");
	}
	size_t bestindeme = 0, worstinpicks = 0;
	double best = deme[0].get_fitness(), worst = pop[picks[0]].get_fitness();
	// Determine the best individual in deme and the worst among the picks.
	for (size_t i = 1; i < Ndeme; ++i) {
		if (deme[i].get_fitness() < best) {
			bestindeme = i;
			best = deme[i].get_fitness();
		}
		if (picks[i] >= pop_size) {
			pagmo_throw(index_error,"pick value exceeds population's size while inserting best in deme");
		}
		if (pop[picks[i]].get_fitness() > worst) {
			worstinpicks = i;
			worst = pop[picks[i]].get_fitness();
		}
	}
	// In place of the worst among the picks, insert the best in deme.
	pop[picks[worstinpicks]] = deme[bestindeme];

}


individual population::checked_individual(const individual &ind) const
{
	try {
		ind.check(*m_problem);
	} catch (const value_error &) {
		const size_t size = ind.get_decision_vector().size();
		if (size != m_problem->getDimension()) {
			pagmo_throw(value_error,"individual's size is incompatible with population");
		} else {
			const individual random_ind = individual(*m_problem);
			std::vector<double> dv(size), v(size);
			for (size_t i = 0; i < size; ++i) {
				if (ind.get_decision_vector()[i] > m_problem->get_ub()[i] ||
				        ind.get_decision_vector()[i] < m_problem->get_lb()[i]) {
					dv[i] = random_ind.get_decision_vector()[i];
					v[i] = random_ind.get_velocity()[i];
				} else {
					dv[i] = ind.get_decision_vector()[i];
					v[i] = ind.get_velocity()[i];
				}
			}
			return individual(*m_problem,dv,v);
		}
	}
	return ind;
}


void population::createRandompopulation(int N)
{
	pop.clear();
	pop.reserve(N);
	for (int i=0; i < N; ++i) {
		pop.push_back(individual(*m_problem));
	}
}


std::ostream &operator<<(std::ostream &s, const population &p)
{
	s << "Problem type: '" << p.problem().id_name() << "'\n";
	for (size_t i = 0; i < p.size(); ++i) {
		s << "Individual #" << i << ": " << p.pop[i].get_fitness() << " " << p.pop[i] << std::endl;
	}
	return s;
}

}