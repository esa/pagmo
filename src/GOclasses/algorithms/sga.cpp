
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

// Created by Dario Izzo on 10/05/08.

#include <cmath>
#include <iostream>
#include <vector>

#include "../../exceptions.h"
#include "../basic/individual.h"
#include "../basic/population.h"
#include "../problems/base.h"
#include "base.h"
#include "sga.h"

// TODO: missing header for rngs.

using namespace std;

namespace pagmo
{
namespace algorithm {

sga::sga(int generationsInit, const double &CRInit, const double &MInit, int insert_bestInit):
		base(),generations(generationsInit),CR(CRInit),M(MInit),insert_best(insert_bestInit),MR(1),MType(2),SType(1),rng(static_rng_uint32()())
{
	if (generationsInit <= 0) {
		pagmo_throw(value_error,"number of generations must be strictly positive");
	}
	if (CRInit < 0) {
		pagmo_throw(value_error,"CR value must be non-negative");
	}
	if (MInit < 0) {
		pagmo_throw(value_error,"M value must be non-negative");
	}
	if (insert_bestInit <= 0) {
		pagmo_throw(value_error,"insert_best value must be positive");
	}
}


sga::sga(int generationsInit, const double &CRInit, const double &MInit, int insert_bestInit, double mutationRange, int mutationType, int selectionType):
		base(),generations(generationsInit),CR(CRInit),M(MInit),insert_best(insert_bestInit),rng(static_rng_uint32()())
{
	if (generationsInit <= 0) {
		pagmo_throw(value_error,"number of generations must be strictly positive");
	}
	if (CRInit < 0) {
		pagmo_throw(value_error,"CR value must be non-negative");
	}
	if (MInit < 0) {
		pagmo_throw(value_error,"M value must be non-negative");
	}
	if (insert_bestInit <= 0) {
		pagmo_throw(value_error,"insert_best value must be positive");
	}
	if (mutationRange <= 0 && mutationType < 2) {
		pagmo_throw(value_error,"mutationRange value must be positive");
	}
	if (mutationType < 0 || mutationType > 2) {
		pagmo_throw(value_error,"mutationType does not corrospond to any defined type");
	}
	if (selectionType < 0 || selectionType > 1) {
		pagmo_throw(value_error,"selectionType does not corrospond to any defined type");
	}
	MR = mutationRange;
	MType = mutationType;
	SType = selectionType;
}

population sga::evolve(const population &deme) const
{
	const problem::base &problem = deme.problem();

	const std::vector<double>& LB = problem.get_lb();
	const std::vector<double>& UB = problem.get_ub();

	int NP = deme.size();
	int D = LB.size();

	vector<double> dummy(D,0);						//used for initialisation purposes
	vector< vector<double> > X(NP,dummy), Xnew(NP,dummy);

	vector<double> fit(NP);							//fitness

	double bestfit=0;
	vector<double> bestX(D,0);

	vector<double> selectionfitness(NP), cumsum(NP), cumsumTemp(NP);
	vector <int> selection(NP);

	int tempID;
	double temp;
	int fitnessID[NP];

	// Initialise the chromosomes and their fitness to that of the initial deme
	for ( int i = 0; i<NP; i++ ) {
		X[i]	=	deme[i].get_decision_vector();
		fit[i]	=	deme[i].get_fitness();
	}

	// Find the best member and store in bestX and bestfit
	bestfit = fit[0];
	bestX = X[0];
	for (int i = 1; i<NP; i++) {		//the int i = 1 jumps the first member as it is already set as the best
		if (fit[i] < bestfit) {
			bestfit = fit[i];
			bestX = X[i];
		}
	}


	// Main SGA loop
	for (int j = 0; j<(int)generations; j++) {
		//for (int k = 0; k < 5; k++)
		//{
		//	cout<<"GENERATIONS: "<<generations<<"  THIS: "<<j<<endl;
		//}

		//1 - Selection
		selectionfitness.clear();
		cumsum.clear();
		cumsumTemp.clear();
		double worstfit=fit[0];
		for (int i = 1; i < NP;i++) {
			if (fit[i] > worstfit) worstfit=fit[i];
		}

		double factor = std::log((double)(j+1))/10;

		// selectionfitness contains a modification of the objective function (something like a fitness) used for the solution ranking
		for (int i = 0; i < NP; i++) {
			selectionfitness.push_back(pow((worstfit - fit[i]),1.0 + factor)); //works for minimisation
		}

		// cumsumTemp contains the cumulative sum of the objective functions in pold
		cumsumTemp.push_back(selectionfitness[0]);
		for (int i = 1; i< NP; i++) {
			cumsumTemp.push_back(cumsumTemp[i - 1] + selectionfitness[i]);
		}

		// cumsum is cumsumTemp normalised so that it increases from 0 to 1.
		for (int i = 0; i < NP; i++) {
			cumsum.push_back(cumsumTemp[i]/cumsumTemp[NP-1]);
		}

		selection.resize(NP,0);

		switch (SType) {
			//0 - select best %20 and reproduce each 5x
		case 0:
			//Sort the individuals according to their fitness
			for (int i=0;i<NP;i++) fitnessID[i]=i;
			for (int i=0; i<(NP-1); ++i) {
				for (int j=i+1; j<NP; ++j) {
					if (selectionfitness[i] < selectionfitness[j]) {
						//swap fitness values
						temp = selectionfitness[i];
						selectionfitness[i] = selectionfitness[j];
						selectionfitness[j] = temp;
						//swap id's
						tempID = fitnessID[i];
						fitnessID[i] = fitnessID[j];
						fitnessID[j] = tempID;
					}
				}
			}
			for (int i = 0;i < NP/5; i++) {
				for (int j = 0; j < 5; j++) {
					selection[(i*5)+j] = fitnessID[i];
				}
			}
			break;

			//1 - roulette selection
		case 1:
			double r2;
			for (int i = 0; i < NP; i++) {
				r2 = drng();
				for (int j = 0; j < NP; j++) {
					if (cumsum[j] > r2) {
						selection[i]=j;
						break;
					}
				}
			}
			break;
		}

		//Xnew is the selected population
		for (int i = 0; i < NP; i++) {
			Xnew[i]=X[selection[i]];
		}

		//2 - Crossover
		{
			int r1,n,L;
			vector<double>  member1,member2;

			for (int i=0; i< NP; i++) {
				//for each member selected from pold (i.e. in pnew)
				member1 = Xnew[i];
				//we select a mating patner (different from the self (i.e. no masturbation))
				do { //FIXME: [MaRu] YOU DON'T DO IT THAT WAY, MAN!!!
					r1 = rng() % NP;
				} while ( r1 == i );
				member2 = Xnew[r1];
				//and we operate crossover
				n = rng() % D;
				L = 0;
				do {
					member1[n] = member2[n];
					n = (n+1) % D;
					L++;
				}  while ( (drng() < CR) && (L < D) );
				Xnew[i] = member1;
			}
		}

		//3 - Mutation
		{
			switch (MType) {

				//Original Gaussian mutation
			case 0:
				for (int i = 0; i < NP;i++) {
					//generate a random mutation vector
					for (int k = 0; k < D;k++)
						if (drng() < M) {
							//Gaussian mutation instead of : Xnew[i][j] = (UB[j]-LB[j])*drng() + LB[j];
							Xnew[i][k] = Xnew[i][k] + sqrt(-2*(std::log(drng())))*cos(2*3.14159265*drng());
							if (Xnew[i][k] > UB[k]) Xnew[i][k] = UB[k];
							if (Xnew[i][k] < LB[k]) Xnew[i][k] = LB[k];
						}
				}
				break;
				//2 - Bounded mutation value
			case 1:
				for (int i = 0; i < NP;i++) {
					for (int j = 0; j < D;j++) {
						if (drng() < M) {
							Xnew[i][j] = (MR-(-MR))*drng()+(-MR);
							if (Xnew[i][j] > UB[j]) Xnew[i][j] = UB[j];
							if (Xnew[i][j] < LB[j]) Xnew[i][j] = LB[j];
						}
					}
				}
				break;
				//3 - Random mutation
			case 2:
				for (int i = 0; i < NP;i++) {
					//generate a random mutation vector
					for (int j = 0; j < D;j++) {
						if (drng() < M) {
							Xnew[i][j] = (UB[j]-LB[j])*drng() + LB[j];
						}
					}
				}
				break;
			}

		}

		//4 - Evaluate Xnew
		{
			for (int i = 0; i < NP;i++) {
				fit[i] = problem.objfun(Xnew[i]);
				if (fit[i] < bestfit) {
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}

		//5 - Reinsert best individual every insert_best generations (TODO: sort the population?)

		if (j % insert_best == 0) {
			int worst=0;
			for (int i = 1; i < NP;i++) {
				if (fit[i] > fit[worst]) worst=i;
			}
			Xnew[worst] = bestX;
			fit[worst] = bestfit;
		}
		X = Xnew;
	} // end of main SGA loop

	//we end by constructing the object population containing the final results
	population popout(problem,0);
	for (int i=0; i<NP; i++) {
		// TODO: WARNING! here we are creating the individuals with empty velocity vectors,
		// we could improve this by, e.g., re-using the velocities from the input population.
		popout.push_back(individual(X[i],std::vector<double>(X[i].size()),fit[i]));
	}
	return popout;
}

void sga::log(std::ostream &s) const
{
	s << "SGA - Generations:" << generations << " CR:" << CR << " M:" << M << " insert_best:" << insert_best <<
	" mutation range:" << MR << " mutation type:" << MType << " selection type:" << SType;
}

}
}
