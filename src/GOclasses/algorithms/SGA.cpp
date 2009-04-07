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

// Created by Dario Izzo on 10/05/08.

#include <cmath>
#include <iostream>
#include <vector>

#include "../../exceptions.h"
#include "../basic/individual.h"
#include "../basic/population.h"
#include "../problems/GOproblem.h"
#include "SGA.h"
#include "go_algorithm.h"

using namespace std;

SGAalgorithm::SGAalgorithm(int generationsInit, const double &CRInit, const double &MInit, int insert_bestInit):
	go_algorithm(),generations(generationsInit),CR(CRInit),M(MInit),insert_best(insert_bestInit),rng(static_rng_uint32()())
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

Population SGAalgorithm::evolve(const Population &deme) const {
	const GOProblem &problem = deme.problem();

    const std::vector<double>& LB = problem.getLB();
    const std::vector<double>& UB = problem.getUB();

    int NP = deme.size();
    int D = LB.size();

	vector<double> dummy(D,0);							//used for initialisation purposes
	vector< vector<double> > X(NP,dummy), Xnew(NP,dummy);

	vector<double> fit(NP);							    //fitness

	double bestfit=0;
	vector<double> bestX(D,0);

	vector<double> selectionfitness(NP), cumsum(NP), cumsumTemp(NP);
	vector <int> selection(NP);


	// Initialise the chromosomes and their fitness to that of the initial deme
	for ( int i = 0; i<NP; i++ ){
			X[i]	=	deme[i].getDecisionVector();
			fit[i]	=	deme[i].getFitness();
   }

   // Find the best member and store in bestX and bestfit
   bestfit = fit[0];
   bestX = X[0];
   for (int i = 1; i<NP; i++){		//the int i = 1 jumps the first member as it is already set as the best
	if (fit[i] < bestfit){
		bestfit = fit[i];
		bestX = X[i];
	}
   }


   // Main SGA loop
   for (int j = 0; j<(int)generations; j++){


		//1 - Selection
			selectionfitness.clear();
			cumsum.clear();
			cumsumTemp.clear();
			double worstfit=fit[0];
			for (int i = 1; i < NP;i++){
				if (fit[i] > worstfit) worstfit=fit[i];
			}

			double factor = std::log((double)(j+1))/10;
			for (int i = 0; i < NP; i++){

				selectionfitness.push_back(pow((worstfit - fit[i]),1.0 + factor)); //works for minimisation

			} // selectionfitness contains a modification of the objective function (something like a fitness) used for the solution ranking

			cumsumTemp.push_back(selectionfitness[0]);
			for (int i = 1; i< NP; i++){
				cumsumTemp.push_back(cumsumTemp[i - 1] + selectionfitness[i]);
			} // cumsumTemp contains the cumulative sum of the objective functions in pold

			for (int i = 0; i < NP; i++){
				cumsum.push_back(cumsumTemp[i]/cumsumTemp[NP-1]);
			} // cumsum is cumsumTemp normalised so that it increases from 0 to 1.

			selection.resize(NP,0);
			double r2;
			for (int i = 0; i < NP; i++){
				r2 = drng();
				for (int j = 0; j < NP; j++){
					if (cumsum[j] > r2){
						selection[i]=j;
						break;
					}
				}
			} //selection contains the roulette selection indexes

			for (int i = 0; i < NP; i++)
			{
				Xnew[i]=X[selection[i]];
			} //Xnew is the selected population



		//2 - Crossover
		{
			int r1,n,L;
			vector<double>  member1,member2;

			for (int i=0; i< NP; i++){
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
			for (int i = 0; i < NP;i++){
				//generate a random mutation vector
				for (int j = 0; j < D;j++)
					if (drng() < M)
					{
						Xnew[i][j] = (UB[j]-LB[j])*drng() + LB[j];
					}
			}
		}

		//4 - Evaluate Xnew
		{
			for (int i = 0; i < NP;i++){
				fit[i] = problem.objfun(Xnew[i]);
				if (fit[i] < bestfit){
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}

		//5 - Reinsert best individual every insert_best generations (TODO: sort the population?)

		if (j % insert_best == 0)
			{
				int worst=0;
				for (int i = 1; i < NP;i++){
					if (fit[i] > fit[worst]) worst=i;
				}
				Xnew[worst] = bestX;
				fit[worst] = bestfit;
			}
		X = Xnew;
   } // end of main SGA loop

   //we end by constructing the object Population containing the final results
   Population popout(problem,0);
   for (int i=0; i<NP; i++){
	popout.push_back(Individual(X[i],std::vector<double>(),fit[i]));
   }
   return popout;
}

void SGAalgorithm::log(std::ostream &s) const
{
	s << "SGA - Generations:" << generations << " CR:" << CR << " M:" << M
		<< " insert_best:" << insert_best;
}
