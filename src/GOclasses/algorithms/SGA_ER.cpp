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
// This algorithm is identical to SGA except for reinsertion. It also evaluated the population in the beginning. This modification has been made to make SGA suitable for no-deterministic/stochastic objective functions, like the ones of Evolutionary Robotics.
#include <cmath>
#include <iostream>
#include <vector>

#include "../../exceptions.h"
#include "../basic/individual.h"
#include "../basic/population.h"
#include "../problems/GOproblem.h"
#include "SGA_ER.h"
#include "go_algorithm.h"

using namespace std;

SGA_ERalgorithm::SGA_ERalgorithm(int generationsInit, const double &CRInit, const double &MInit, int insert_bestInit, int initial_seed):
	go_algorithm(),generations(generationsInit),CR(CRInit),M(MInit),insert_best(insert_bestInit),rng(uint32_t(initial_seed))
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
	drng.seed(uint32_t(initial_seed));
}

Population SGA_ERalgorithm::evolve(const Population &deme_in) const {
	Population deme = deme_in;

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


//Evaluate the population again as the random seed has changed
	for (int i = 0; i < NP;i++){
		deme.set_individual(i,Individual(problem, deme[i].getDecisionVector()));
 	}

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


   // Main SGA_ER loop
for (int j = 0; j<(int)generations; j++){   


		//1 - Selection
			selectionfitness.clear();
			cumsum.clear();
			cumsumTemp.clear();
			double worstfit=fit[0];
			for (int i = 1; i < NP;i++){
				if (fit[i] > worstfit) worstfit=fit[i];
			}

			//double factor = std::log((double)(j+1))/10;
			for (int i = 0; i < NP; i++){

				selectionfitness.push_back((worstfit - fit[i])); //works for minimisation
				//WE DONT NEED TO LOG THE SELECTION IN ER EXPERIMENTS!?!
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
				//BUG? WHEREVER THERE WAS j IN THE GENERATION LOOP I REPLACED IT WITH k
				for (int k = 0; k < NP; k++){
					if (cumsum[k] > r2){
						selection[i]=k;
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
				for (int k = 0; k < D;k++)
					if (drng() < M)
					{
					  //Gaussian mutation instead of : Xnew[i][j] = (UB[j]-LB[j])*drng() + LB[j];
					  Xnew[i][k] = Xnew[i][k] + sqrt(-2*(std::log(drng())))*cos(2*3.14159265*drng());
					  if (Xnew[i][k] > UB[k])
					    Xnew[i][k] = UB[k];
					 if (Xnew[i][k] < LB[k])
					    Xnew[i][k] = LB[k];
					}
			}
		}

// 		//4 - Evaluate Xnew
// 		{
		if (generations>1&&j!=0){
		  for (int i = 0; i < NP;i++){
 				fit[i] = problem.objfun(Xnew[i]);
			}
			bestfit = fit[0];
			bestX = X[0];
   for (int i = 1; i<NP; i++){		//the int i = 1 jumps the first member as it is already set as the best
	if (fit[i] < bestfit){
		bestfit = fit[i];
		bestX = X[i];
	}
   }
		}
// 				if (fit[i] < bestfit){
// 					bestfit = fit[i];
// 					bestX = Xnew[i];
// 				}
// 			}
// 		}

		//5 - Elitism every insert_best generations - Reinsert best individual every insert_best generations (TODO: sort the population?)
		//-always reinsert the best - even if you have a bunch of generatios running together
		//if ((j+1) % insert_best == 0)
		  //{
				int worst=0;
				for (int i = 1; i < NP;i++){
					if (fit[i] > fit[worst]) worst=i;
				}
				Xnew[worst] = bestX;
				fit[worst] = bestfit;
				//}
		X = Xnew;
	 } // end of main SGA_ER loop

   
//we end by constructing the object Population containing the final results
   Population popout(problem,0);
   for (int i=0; i<NP; i++){
	popout.push_back(Individual(X[i],std::vector<double>(),fit[i]));
   }
   return popout;

}

void SGA_ERalgorithm::log(std::ostream &s) const
{
	s << "SGA_ER - Generations:" << generations << " CR:" << CR << " M:" << M
		<< " insert_best:" << insert_best;
}
