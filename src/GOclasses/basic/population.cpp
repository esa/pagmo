/*
 *  population.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include <cmath>
#include <functional>

#include "GOproblem.h"
#include "population.h"
#include "rng.h"

	Population::Population(GOProblem &p, int N) {
		createRandomPopulation(p,N);
	}

	void Population::createRandomPopulation(GOProblem &problem, int N){
		pop.clear();
		for (int i=0; i < N; i++){
			pop.push_back(Individual(problem));
		}
	};

	void Population::push_back(const Individual &x){
		pop.push_back(x);
	};

	size_t Population::size() const {
		return pop.size();
	};

	const Individual &Population::extractBestIndividual() const {
		return extract_most<std::less<double> >();
	}

	const Individual &Population::extractWorstIndividual() const {
		return extract_most<std::greater<double> >();
	}

	Population Population::extractRandomDeme(int N, std::vector<int> &picks, rng_double &drng){
		Population deme;
		std::vector<int> PossiblePicks;
		int Pick;

		for (unsigned int i=0; i < pop.size(); i++)
			PossiblePicks.push_back(i);


		for (int i=0; i < N; i++){
			//we pick a random position between 0 and popsize-1
			Pick = (int)(drng() * PossiblePicks.size());
			//and store it
			picks.push_back(PossiblePicks[Pick]);
			//we insert the corresponding individual in the deme
			deme.push_back(pop[PossiblePicks[Pick]]);
			//and erase it from the possible picks
			PossiblePicks.erase(PossiblePicks.begin() + Pick);
		}
		return deme;
	};

	void Population::insertDeme(const Population &deme, const std::vector<int> &picks){
		for (unsigned int i=0; i<picks.size(); i++){
			if ( deme[i].getFitness() < pop[picks[i]].getFitness() ){
				pop[picks[i]] = deme[i];
			}
		}
	}

	void Population::insertDemeForced(const Population &deme, const std::vector<int> &picks){
		for (unsigned int i=0; i<picks.size(); i++){
				pop[picks[i]] = deme[i];
		}
	}

	void Population::insertBestInDeme(const Population &deme, const std::vector<int> &picks){
		const int Ndeme = deme.size();

		int bestindeme =  0;
		int worstinpicks = 0;
		double best = deme[0].getFitness();
		double worst = pop[picks[0]].getFitness();

		for (int i=1; i<Ndeme; i++){
			if ( deme[i].getFitness() < best){
				bestindeme = i;
				best = deme[i].getFitness();
			}
			if ( pop[picks[i]].getFitness() > worst) {
				worstinpicks = i;
				worst = pop[picks[i]].getFitness();
			}
		}
		pop[picks[worstinpicks]] = deme[bestindeme];

	}

	double Population::evaluateMean() const {
		double mean=0;
		int size = 0;
		size = pop.size();
		for (int i=0; i<size; i++){
			mean=mean+pop[i].getFitness();
		}
		mean = mean/(double)size;
		return mean;
	}

	double Population::evaluateStd() const {
		double Std=0,mean=0;
		int size = 0;
		size = pop.size();
		mean = evaluateMean();

		for (int i=0; i<size; i++){
			Std=Std+pow((pop[i].getFitness()-mean),2.0);
		}
		Std = sqrt(Std/size);
		return Std;
	}

	Individual &Population::operator[](const size_t &index){
		return pop[index];
	};

	const Individual &Population::operator[](const size_t &index) const{
		return pop[index];
	};

	std::ostream& operator<<(std::ostream& s, Population& pop){
	for (unsigned int i=0;i<pop.size(); i++){
		s << "Individual #" << i << ": " << pop[i].getFitness() << " " << pop[i] << std::endl;
	}
	return s;
    }

