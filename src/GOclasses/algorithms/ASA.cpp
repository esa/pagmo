/*
 *  ASA.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include <boost/thread/thread.hpp>
#include <cmath>
#include <iostream>

#include "../../exceptions.h"
#include "../basic/individual.h"
#include "../basic/population.h"
#include "ASA.h"
#include "go_algorithm.h"

	ASAalgorithm::ASAalgorithm(int niterTotInit, const GOProblem &problem, const double &Ts, const double &Tf):
		go_algorithm(problem) {
		if (niterTotInit < 0) {
			pagmo_throw(value_error,"number of generations must be nonnegative");
		}
		if (Ts <= 0 || Tf <= 0 || Ts <= Tf) {
			pagmo_throw(value_error,"temperatures must be positive and Ts must be greater than Tf");
		}
		niterTot = niterTotInit;
		niterTemp = 1;
		niterRange = 20;
		niterOuter = niterTot / (niterTemp * niterRange * SolDim);
		T0 = Ts;
		Tcoeff = std::pow(Tf/Ts,1.0/(double)(niterOuter));
		StartStep = 1;
	}

Population ASAalgorithm::evolve(const Population &pop, GOProblem &problem) {
	if (pop.size() == 0) {
		pagmo_throw(index_error,"population's size cannot be null");
	}

    const std::vector<double> &LB = problem.getLB();
    const std::vector<double> &UB = problem.getUB();
    const Individual &x0 = pop[0];

	std::vector<double> xNEW = x0.getDecisionVector(), xOLD = xNEW;
	double fNEW = x0.getFitness(), fOLD = fNEW;
	if (xNEW.size() != SolDim) {
		pagmo_throw(value_error,"discrepancy between individual size and problem size.");
	}
	std::vector<double> Step(SolDim,StartStep);
	std::vector<int> acp(SolDim,0) ;
	double ratio = 0, currentT = T0, prob = 0,  r = 0;

	//Main SA loops
	for (size_t jter = 0; jter < niterOuter; ++jter) {
		for (size_t mter = 0; mter < niterTemp; ++mter) {
			for (size_t kter = 0; kter < niterRange; ++kter) {
				size_t nter =  (size_t)(drng() * SolDim);
				for (size_t numb = 0; numb < SolDim ; ++numb) {
					nter = (nter + 1) % SolDim;
					//We modify the current point actsol by mutating its nter component within
					//a Step that we will later adapt
					r = 2.0 * drng() - 1.0; //random number in [-1,1]
					xNEW[nter] = xOLD[nter] + r * Step[nter] * ( UB[nter] - LB[nter] );

					// If new solution produced is infeasible ignore it
					if ((xNEW[nter] > UB[nter]) || (xNEW[nter] < LB[nter])) {
						xNEW[nter]=xOLD[nter];
						continue;
					}
					//And we valuate the objective function for the new point
					fNEW = problem.objfun(xNEW);

					// We decide wether to accept or discard the point
					if (fNEW < fOLD) {
						//accept
						xOLD[nter] = xNEW[nter];
						fOLD = fNEW;
						acp[nter]++;	//Increase the number of accepted values
					} else {
						//test it with Boltzmann to decide the acceptance
						prob = exp ( (fOLD - fNEW )/ currentT );

						// we compare prob with a random probability.
						if (prob > drng()) {
							xOLD[nter] = xNEW[nter];
							fOLD = fNEW;
							acp[nter]++;	//Increase the number of accepted values
						} else {
							xNEW[nter] = xOLD[nter];
						}
					} // end if
				} // end for(nter = 0; ...
			} // end for(kter = 0; ...


			// adjust the step (adaptively)
			for (size_t iter = 0; iter < SolDim; ++iter) {
				ratio = (double)acp[iter]/(double)niterRange;
				acp[iter] = 0;  //reset the counter
				if (ratio > .6) {
					//too many acceptances, increase the step by a factor 3 maximum
					Step[iter] = Step [iter] * (1 + 2 *(ratio - .6)/.4);
				} else {
					if (ratio < .4) {
						//too few acceptance, decrease the step by a factor 3 maximum
						Step [iter]= Step [iter] / (1 + 2 * ((.4 - ratio)/.4));
					};
				};
				//And if it becomes too large, reset it to its initial value
				if ( Step[iter] > StartStep ) {
						Step [iter] = StartStep;
				};
			}
		}
		// Cooling schedule
		currentT *= Tcoeff;
	}

	Population newpop;
	if (fOLD < x0.getFitness()){
		newpop.push_back(Individual(xOLD,x0.getVelocity(),fOLD));
	} else {
		newpop.push_back(x0);
	}
	return newpop;
	}

void mt_asa_algorithm::evolve(const Population &pop, GOProblem &prob) {
	boost::thread thr(evolver(),m_asa,pop,&prob);
}
