/*
 *  DE.cpp
 *  SeGMO
 *
 *  Created by Dario Izzo on 5/18/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 * Algorithm adapted from de36.c available to download
 * from http://www.icsi.berkeley.edu/~storn/code.html#csou
 */

#include <iostream>
#include <vector>

#include "DE.h"

using namespace std;

void DEalgorithm::initDE(int generationsInit, int SolDimInit, double FInit, double CRInit, int strategyInit, uint32_t randomSeed){
	generations = generationsInit;
	strategy = strategyInit;
	F = FInit;
	CR = CRInit;
	SolDim = SolDimInit;
	drng.seed(randomSeed);
}

Population DEalgorithm::evolve(Population deme, GOProblem& problem){

    const std::vector<double>& LB = problem.getLB();
    const std::vector<double>& UB = problem.getUB();

	int NP = deme.size();
	int D = SolDim;

	vector<double> dummy(D),tmp(D);						//dummy is used for initialisation purposes, tmp to contain
														//the mutated candidate
	vector< vector<double> > popold(NP,dummy),popnew(NP,dummy), popswap(NP,dummy);
	vector<double> fit(NP);								//chromosome fitness

	double newfitness;									//new fitness of the mutaded candidate
	double gbfit;										//global best fitness
	vector<double> gbX(D);								//global best chromosome
	vector<double> gbIter(D);							//best chromosome of current iteration

	int i,gen,L,n;										//counters


	// Initialise the chromosome (class Individual) values, and their fitness to that of the deme
	for (i = 0; i<NP; i++){
			popold[i] = deme[i].getDecisionVector();
			popnew[i] = deme[i].getDecisionVector();
			fit[i] = deme[i].getFitness();
   }

   // Initialise the global bests
   gbX=popold[0];
   gbfit=fit[0];

   for (i = 1; i<NP; i++){		//the int i = 1 jumps the first member as it is already set as the best
	if (fit[i] < gbfit){
		gbfit = fit[i];			// save best member ever
		gbX = popold[i];
	}
   }
   gbIter = gbX;				// save best member of generation

   // Main DE iterations
   int r1,r2,r3,r4,r5;	//indexes to the selected population members
   for (gen=0; gen < generations; gen++){
		//Start of the loop through the deme
		for (i=0; i<NP; i++){
			do{                        /* Pick a random population member */
			                          /* Endless loop for NP < 2 !!!     */
				r1 = (int)(drng()*NP);
			}while(r1==i);

			do{                        /* Pick a random population member */
									   /* Endless loop for NP < 3 !!!     */
				r2 = (int)(drng()*NP);
			}while((r2==i) || (r2==r1));

			do{                        /* Pick a random population member */
									   /* Endless loop for NP < 4 !!!     */
				r3 = (int)(drng()*NP);
			}while((r3==i) || (r3==r1) || (r3==r2));

			do{                        /* Pick a random population member */
									   /* Endless loop for NP < 5 !!!     */
				r4 = (int)(drng()*NP);
			}while((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));

			do{                        /* Pick a random population member */
									   /* Endless loop for NP < 6 !!!     */
				r5 = (int)(drng()*NP);
			}while((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));


/*-------DE/best/1/exp--------------------------------------------------------------------*/
/*-------Our oldest strategy but still not bad. However, we have found several------------*/
/*-------optimization problems where misconvergence occurs.-------------------------------*/
			if (strategy == 1){ /* strategy DE0 (not in our paper) */
				tmp = popold[i];
				n = (int)(drng()*D);
				L = 0;
				do{
					tmp[n] = gbIter[n] + F*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%D;
					L++;
				}while((drng() < CR) && (L < D));
			}

/*-------DE/rand/1/exp-------------------------------------------------------------------*/
/*-------This is one of my favourite strategies. It works especially well when the-------*/
/*-------"gbIter[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5---------*/
/*-------as a first guess.---------------------------------------------------------------*/
			else if (strategy == 2){ /* strategy DE1 in the techreport */
				tmp = popold[i];
				n = (int)(drng()*D);
				L = 0;
				do{
					tmp[n] = popold[r1][n] + F*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%D;
					L++;
				}while((drng() < CR) && (L < D));
			}

/*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
/*-------This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.------*/
/*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
/*-------should play around with all three control variables.----------------------------*/
			else if (strategy == 3){ /* similiar to DE2 but generally better */
				tmp = popold[i];
				n = (int)(drng()*D);
				L = 0;
				do{
					tmp[n] = tmp[n] + F*(gbIter[n] - tmp[n]) + F*(popold[r1][n]-popold[r2][n]);
					n = (n+1)%D;
					L++;
				}while((drng() < CR) && (L < D));
			}
/*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
			else if (strategy == 4){
				tmp = popold[i];
				n = (int)(drng()*D);
				L = 0;
				do{
					tmp[n] = gbIter[n] +
					(popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					n = (n+1)%D;
					L++;
				}while((drng() < CR) && (L < D));
			}
/*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
			else if (strategy == 5){
				tmp = popold[i];
				n = (int)(drng()*D);
				L = 0;
				do{
					tmp[n] = popold[r5][n] +
					(popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					n = (n+1)%D;
					L++;
				}while((drng() < CR) && (L < D));
			}

/*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

/*-------DE/best/1/bin--------------------------------------------------------------------*/
			else if (strategy == 6){
				tmp = popold[i];
				n = (int)(drng()*D);
					for (L=0; L<D; L++){ /* perform D binomial trials */
						if ((drng() < CR) || L == (D-1)){ /* change at least one parameter */
							tmp[n] = gbIter[n] + F*(popold[r2][n]-popold[r3][n]);
						}
					n = (n+1)%D;
					}
			}
/*-------DE/rand/1/bin-------------------------------------------------------------------*/
			else if (strategy == 7){
				tmp = popold[i];
				n = (int)(drng()*D);
				for (L=0; L<D; L++){ /* perform D binomial trials */
					if ((drng() < CR) || L == (D-1)){ /* change at least one parameter */
						tmp[n] = popold[r1][n] + F*(popold[r2][n]-popold[r3][n]);
					}
					n = (n+1)%D;
				}
			}
/*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
			else if (strategy == 8){
				tmp = popold[i];
				n = (int)(drng()*D);
				for (L=0; L<D; L++){ /* perform D binomial trials */
					if ((drng() < CR) || L == (D-1)){ /* change at least one parameter */
						tmp[n] = tmp[n] + F*(gbIter[n] - tmp[n]) + F*(popold[r1][n]-popold[r2][n]);
					}
					n = (n+1)%D;
				}
			}
/*-------DE/best/2/bin--------------------------------------------------------------------*/
			else if (strategy == 9){
				tmp = popold[i];
				n = (int)(drng()*D);
				for (L=0; L<D; L++){ /* perform D binomial trials */
					if ((drng() < CR) || L == (D-1)){ /* change at least one parameter */
						tmp[n] = gbIter[n] +
						(popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					}
					n = (n+1)%D;
				}
			}
/*-------DE/rand/2/bin--------------------------------------------------------------------*/
			else if (strategy == 10){
				tmp = popold[i];
				n = (int)(drng()*D);
				for (L=0; L<D; L++){ /* perform D binomial trials */
					if ((drng() < CR) || L == (D-1)){ /* change at least one parameter */
						tmp[n] = popold[r5][n] +
						(popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					}
					n = (n+1)%D;
				}
			}


/*=======Trial mutation now in tmp[]. Force feasibility and how good this choice really was.==================*/
			// a) feasibility
			int i2;
			i2=0;
			while  (i2<D) {
				if ((tmp[i2] < LB[i2]) || (tmp[i2] > UB[i2]))
					tmp[i2] = drng()*(UB[i2]-LB[i2]) + LB[i2];
				i2++;
			}
			newfitness = problem.objfun(tmp);    /* Evaluate new vector in tmp[] */
			//cout << "\n" << newfitness << endl;
			//for (int g=0;g<NP;g++){
			//cout << tmp[g] << " ";
			//}
			//b) how good?
			if (newfitness <= fit[i]){   /* improved objective function value ? */
				fit[i]=newfitness;
				popnew[i] = tmp;
				if (newfitness<gbfit){         /* Was this a new minimum for the deme? */
											   /* if so...*/
					gbfit=newfitness;          /* reset gbfit to new low...*/
					gbX=tmp;
				}
			}
			else {
				popnew[i] = popold[i];
			}
			/* swap population arrays. New generation becomes old one */

		}//End of the loop through the deme

		/* Save best population member of current iteration */
		gbIter = gbX;

		/* swap population arrays. New generation becomes old one */
		for (int ii=0;ii<NP;ii++){
			popswap[ii] = popold[ii];
			popold[ii] = popnew[ii];
			popnew[ii] = popswap[ii];
		}

   }//end main DE iterations

   //we end by constructing the object Population containing the final results
   Population popout;
   Individual dummy2;
   vector <double> Xini(D),Vfin(D);
   for (i=0; i<NP; i++){
	dummy2.setDecisionVector(popold[i]);
	dummy2.setFitness(fit[i]);
	Xini = deme[i].getDecisionVector();
	for (int j=0; j<D; j++){
			Vfin[j] = popold[i][j] - Xini[j];
	}
	//dummy2.setVelocity(Vfin);
	dummy2.setVelocity(deme[i].getVelocity());
	//X[i] - deme[i].getDecisionVector());  DOES NOT WORK AS VECTOR CLASS DOES NOT ACCEPT MINUS AS OPERATOR

	popout.addIndividual(dummy2);
   }
   return popout;
}






