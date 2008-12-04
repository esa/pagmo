/*
 *  individual.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 Àdvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include "individual.h"
#include "PkRandom.h"

	void Individual::createRandomIndividual(std::vector<double> LB, std::vector<double> UB, Pk::Random32& rng){

	    //We first delete the vector content if any
		x.clear();
		v.clear();

		//We then push back random numbers to fill a new random chromosome
		double dummy;
		for (unsigned int i=0; i < LB.size(); i++){
			dummy = LB[i] + Pk::nextDouble(rng) * (UB[i] - LB[i]);
			x.push_back(dummy);
		}
		//And a random velocity
		for (unsigned int i=0; i < LB.size(); i++){
			dummy = Pk::nextDouble(rng) * (UB[i] - LB[i]);  //initial velocity range
			v.push_back(dummy);
		}
	};//createRandomIndividual

	void Individual::resetVelocity(std::vector<double> LB, std::vector<double> UB, Pk::Random32& rng){
		v.clear();
		double dummy;
		for (unsigned int i=0; i < LB.size(); i++){
			dummy = (2*Pk::nextDouble(rng) - 1) * (UB[i] - LB[i]);
			v.push_back(dummy);
		}
	}

	double Individual::evaluateFitness(GOProblem& problem){
		this->fitness = problem.objfun(x);
		return this->fitness;
	};

	double Individual::getFitness() const {
		return fitness;
	};

	void Individual::setFitness(double fitnessnew){
		fitness = fitnessnew;
	};

	std::vector<double> Individual::getDecisionVector() const{
		return x;
	};

	void Individual::setDecisionVector(std::vector<double> xnew){
		x = xnew;
	};

	std::vector<double> Individual::getVelocity() const {
		return v;
	};

	void Individual::setVelocity(std::vector<double> vnew){
		v = vnew;
	};

	double& Individual::operator[](int index) {
		return x[index];
	};

	void Individual::operator=(const Individual& newindividual){
			x		=	newindividual.getDecisionVector();
			v		=	newindividual.getVelocity();
			fitness =	newindividual.getFitness();
	};

    std::ostream& operator<<(std::ostream& s, Individual& x){
        for (unsigned int i=0;i<x.getDecisionVector().size(); i++){
            s << x[i] << " ";
        }
        return s;
    }

