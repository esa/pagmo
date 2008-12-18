/*
 *  ASA.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include "ASA.h"
#include "PkRandom.h"
#include <iostream>

Population ASAalgorithm::evolve(Individual x0, GOProblem& problem) {

    const std::vector<double>& LB = problem.getLB();
    const std::vector<double>& UB = problem.getUB();

	std::vector<double> xNEW = x0.getDecisionVector();
	std::vector<double> xOLD = xNEW;
	double fNEW;
	fNEW = x0.getFitness();
	double fOLD;
	fOLD = fNEW;
	SolDim = xNEW.size();
	std::vector<double> Step(SolDim,StartStep);
	std::vector<int> acp(SolDim,0) ;
	double ratio=0;
	double currentT = T0;
	double prob=0;
	double r=0;

	//Main SA loops
	for ( int jter=0; jter < niterOuter; jter++){
		for ( int mter = 0; mter < niterTemp; mter++){
			for ( int kter = 0 ; kter < niterRange; kter++){
				int nter =  (int)(Pk::nextDouble(rng)*SolDim);
				for ( int numb = 0; numb < SolDim ;numb++){
					nter=(nter+1) % SolDim;
					//We modify the current point actsol by mutating its nter component within
					//a Step that we will later adapt
					r = 2.0*Pk::nextDouble(rng)-1.0; //random number in [-1,1]
					xNEW[nter] = xOLD[nter] + r * Step[nter] * ( UB[nter] - LB[nter] );

					// If new solution produced is infeasible ignore it
					if ((xNEW[nter] > UB[nter]) || (xNEW[nter] < LB[nter])){
						xNEW[nter]=xOLD[nter];
						continue;
					}
					//And we valuate the objective function for the new point
					fNEW = problem.objfun(xNEW);

					// We decide wether to accept or discard the point
					if (fNEW < fOLD)  //accept
					{
						xOLD[nter] = xNEW[nter];
						fOLD = fNEW;
						acp[nter]++;	//Increase the number of accepted values
					}
					else  //test it with Boltzmann to decide the acceptance
					{

						prob = exp ( (fOLD - fNEW )/ currentT );

						// we compare prob with a random probability.
						if (prob > Pk::nextDouble(rng))
						{
							xOLD[nter] = xNEW[nter];
							fOLD = fNEW;
							acp[nter]++;	//Increase the number of accepted values
						}
						else{
							xNEW[nter] = xOLD[nter];
						}
					} // end if
				} // end for(nter = 0; ...
			} // end for(kter = 0; ...


			// adjust the step (adaptively)
			for (int iter=0; iter < SolDim; iter++)
			{
				ratio = (double)acp[iter]/(double)niterRange;
				acp[iter]=0;  //reset the counter
				if (ratio > .6) //too many acceptances, increase the step by a factor 3 maximum
				{
					Step[iter] = Step [iter] * (1 + 2 *(ratio - .6)/.4);
				}
				else
				{
					if (ratio < .4) //too few acceptance, decrease the step by a factor 3 maximum
					{
						Step [iter]= Step [iter] / (1 + 2 * ((.4 - ratio)/.4));
					};
				};
				//And if it becomes too large, reset it to its initial value
				if ( Step[iter] > StartStep ){
						Step [iter] = StartStep;
				};
			}
		}
		// Cooling schedule
		currentT *= Tcoeff;
	}

	Population newpop;
	if (fOLD < x0.getFitness()){
	x0.setDecisionVector( xOLD );
	x0.setFitness( fOLD );
	}
	newpop.addIndividual( x0 );
	return newpop;
	}

	void ASAalgorithm::initASA(int niterTotInit, int niterTempInit, int niterRangeInit, int SolDimInit, double T0Init, double TcoeffInit, double StartStepInit, uint32_t randomSeed){
		niterTot=niterTotInit;
		niterTemp=niterTempInit;
		niterRange=niterRangeInit;
		SolDim=SolDimInit;
		T0=T0Init;
		Tcoeff=TcoeffInit;
		StartStep=StartStepInit;
		niterOuter = niterTot / (niterTemp * niterRange * SolDim);
		rng.seed(randomSeed);
	}

	void ASAalgorithm::initASA(int niterTotInit,int SolDimInit, double Ts, double Tf, uint32_t randomSeed){
		niterTot=niterTotInit;
		niterTemp=1;
		niterRange=20;
		SolDim=SolDimInit;
		niterOuter = niterTot / (niterTemp * niterRange * SolDim);
		T0=Ts;
		Tcoeff=pow(Tf/Ts,1.0/(double)(niterOuter));
		StartStep=1;
		rng.seed(randomSeed);
	}
