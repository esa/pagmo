/*
 *  SGA.cpp
 *  PaGMO
 *
 *  Created by Dario Izzo on 10/5/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>

#include "SGA.h"

using namespace std;

void SGAalgorithm::initSGA(int generationsInit, int SolDimInit, double CRInit, double MInit, int insert_bestInit, uint32_t randomSeed){
	generations = generationsInit;
	SolDim = SolDimInit;
	CR = CRInit;
	M = MInit;
	insert_best = insert_bestInit;
	rng.seed(randomSeed);
	drng.seed(randomSeed);
}

Population SGAalgorithm::evolve(Population deme, GOProblem& problem){

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
   for (int j = 0; j<generations; j++){


		//1 - Selection
			selectionfitness.clear();
			cumsum.clear();
			cumsumTemp.clear();
			double worstfit=fit[0];
			for (int i = 1; i < NP;i++){
				if (fit[i] > worstfit) worstfit=fit[i];
			}

			double factor = log((double)(j+1))/10;
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
   Population popout;
   Individual dummy2;
   for (int i=0; i<NP; i++){
	dummy2.setDecisionVector(X[i]);
	dummy2.setFitness(fit[i]);
	//dummy2.setVelocity(V[i]);
	popout.addIndividual(dummy2);
   }
   return popout;
}




