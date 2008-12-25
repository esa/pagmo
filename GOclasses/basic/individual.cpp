/*
 *  individual.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 Àdvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include "individual.h"
#include "rng.h"

	Individual::Individual(const std::vector<double> &LB, const std::vector<double> &UB, rng_double_type &drng):
		x(LB.size()),v(LB.size()),fitness(0) {
		//We then push back random numbers to fill a new random chromosome
		const size_t size = LB.size();
		for (size_t i = 0; i < size; ++i){
			x[i] = LB[i] + drng() * (UB[i] - LB[i]);
			v[i] = drng() * (UB[i] - LB[i]);
		}
	}

	void Individual::resetVelocity(const std::vector<double> &LB, const std::vector<double> &UB, rng_double_type &drng){
		v.clear();
		double dummy;
		for (unsigned int i=0; i < LB.size(); i++){
			dummy = (2*drng() - 1) * (UB[i] - LB[i]);
			v.push_back(dummy);
		}
	}

	double Individual::evaluateFitness(GOProblem &problem){
		this->fitness = problem.objfun(x);
		return this->fitness;
	};

	double Individual::getFitness() const {
		return fitness;
	};

	const std::vector<double> &Individual::getDecisionVector() const{
		return x;
	};

	void Individual::setDecisionVector(const std::vector<double> &xnew){
		x = xnew;
	};

	const std::vector<double> &Individual::getVelocity() const {
		return v;
	};

	void Individual::setVelocity(const std::vector<double> &vnew){
		v = vnew;
	};

	double &Individual::operator[](int index) {
		return x[index];
	};

	const double &Individual::operator[](int index) const {
		return x[index];
	};

	void Individual::operator=(const Individual &newindividual){
			x		=	newindividual.getDecisionVector();
			v		=	newindividual.getVelocity();
			fitness =	newindividual.getFitness();
	};

    std::ostream& operator<<(std::ostream &s, const Individual &x){
        for (unsigned int i=0;i<x.getDecisionVector().size(); i++){
            s << x[i] << " ";
        }
        return s;
    }

