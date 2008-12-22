/*
 *  PSO.cpp
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 Advanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#include <vector>

#include "PSO.h"

using namespace std;

void PSOalgorithm::initPSO(int generationsInit, int SolDimInit, double omegaInit, double eta1Init, double eta2Init,double vcoeffInit, uint32_t randomSeed){
	generations = generationsInit;
	SolDim = SolDimInit;
	omega = omegaInit;
	eta1 = eta1Init;
	eta2 = eta2Init;
	vcoeff = vcoeffInit;
	drng.seed(randomSeed);
}

Population PSOalgorithm::evolve(Population deme, GOProblem& problem){

    const std::vector<double>& LB = problem.getLB();
    const std::vector<double>& UB = problem.getUB();

	int n = deme.size();
	int m = LB.size();

	vector<double> dummy(m,0);							//used for initialisation purposes
	vector< vector<double> > X(n,dummy);
	vector< vector<double> > V(n,dummy);

	vector<double> fit(n);							//particle fitness

	double gbfit;										//global best fitness
	vector<double> gbX(m);								//global best chromosome

	vector<double> lbfit(n);							//local best fitness
	vector< vector<double> > lbX(n,dummy);				//local best chromosome

	double vwidth;										//Width of the search space
    vector<double> MINV(m),MAXV(m);						//Maximum and minumum velocity allowed

	// Initialise the particle (class Individual) positions, their velocities and their fitness to that of the deme
	for ( int i = 0; i<n; i++ ){
			X[i]	=	deme[i].getDecisionVector();
			V[i]	=	deme[i].getVelocity();
			fit[i]	=	deme[i].getFitness();
   }

   // Initialise the minimum and maximum velocity
   for ( int i = 0; i<m; i++ ) {
		vwidth = (UB[i]-LB[i])/vcoeff;
		MINV[i] = -1.0*vwidth;
		MAXV[i] = vwidth;
	}

   // Initialise the global and local bests
   gbX=X[0];
   gbfit=fit[0];

   lbX=X;			//at the first generation the local best position is the particle position
   lbfit=fit;		//same for the fitness

   for (int i = 1; i<n; i++){		//the int i = 1 jumps the first member as it is already set as the best
	if (fit[i] < gbfit){
		gbfit = fit[i];
		gbX = X[i];
	}
   }


   // Main PSO loop
   for (int j = 0; j<generations; j++){


		//1 - move the particles and check that velocity and positions are in allowed range
		for (int ii = 0; ii< n; ii++){
			for (int jj = 0; jj< m; jj++){

				//new velocity
				V[ii][jj] = omega * V[ii][jj] + eta1 * drng() * (lbX[ii][jj] - X[ii][jj]) + eta2 * drng() * (gbX[jj] - X[ii][jj]);

				//check that it is within the allowed velocity range
				if ( V[ii][jj] > MAXV[jj] )
					V[ii][jj] = MAXV[jj];

				else if ( V[ii][jj] < MINV[jj] )
					V[ii][jj] = MINV[jj];

				//new position
				X[ii][jj] = X[ii][jj] + V[ii][jj];

				if (X[ii][jj] < LB[jj])
					X[ii][jj] = drng() * (UB[jj] - LB[jj]) + LB[jj];
				else if (X[ii][jj] > UB[jj])
					X[ii][jj] = drng() * (UB[jj] - LB[jj]) + LB[jj];
			}

			//We evaluate the new individual fitness now as to be able to update immediately the global best
			fit[ii] = problem.objfun(X[ii]);
			//update local and global best
			if (fit[ii] < lbfit[ii]){
				lbfit[ii] = fit[ii];	//local best
				lbX[ii] = X[ii];
				if(fit[ii] < gbfit){
					gbfit = fit[ii];	//global best
					gbX	= X[ii];
					ii=ii;
				}
			}
		} //End of loop on the population members


   } // end of main PSO loop

   //we end by constructing the object Population containing the final results
   Population popout;
   Individual dummy2;
   for (int i=0; i<n; i++){
	dummy2.setDecisionVector(lbX[i]);
	dummy2.setFitness(lbfit[i]);
	dummy2.setVelocity(V[i]);
	popout.addIndividual(dummy2);
   }
   return popout;
}






