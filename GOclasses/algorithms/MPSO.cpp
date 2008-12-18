/*
 *  MPSO.cpp
 *  PaGMO
 *
 *  Created by Dario Izzo on 10/23/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "MPSO.h"
#include "vector"

using namespace std;

void MPSOalgorithm::initMPSO(int generationsInit, int SolDimInit, double omegaInit, double eta1Init, double eta2Init,double vcoeffInit, int nswarmsInit, uint32_t randomSeed){
	generations = generationsInit;
	SolDim = SolDimInit;
	omega = omegaInit;
	eta1 = eta1Init;
	eta2 = eta2Init;
	vcoeff = vcoeffInit;
	nswarms = nswarmsInit;
	rng.seed(randomSeed);
}

Population MPSOalgorithm::evolve(Population deme, GOProblem& problem){

    const std::vector<double>& LB = problem.getLB();
    const std::vector<double>& UB = problem.getUB();

	int NP = deme.size()/nswarms;							//potentially dangerous when the deme size is not divisible by the numberof swamrs
	int D = LB.size();

	vector<double> dummy(D,0);								//used for initialisation purposes
	vector<double> dummy2(NP,0);							//used for initialisation purposes
	vector< vector<double> > dummy3(NP,dummy);				//used for initialisation purposes

	//Particle position: X[i][j][k] = i-th swarm, j-th individual, k-th phenotype
	vector< vector< vector<double> > > X(nswarms,dummy3);
	//Particle velocity: V[i][j][k] = i-th swarm, j-th individual, k-th phenotype
	vector< vector< vector<double> > > V(nswarms,dummy3);
	//Particle fitness: fit[i][j] = i-th swarm, j-th individual
	vector< vector<double> > fit(nswarms,dummy2);
	//Global best swarm fitness: gbfit[i] = i-th swarm
	vector<double> gbfit(nswarms,0);
	//Global best swarm individual: gbX[i][j] = i-th swarm, j-th phenotype
	vector< vector<double> > gbX(nswarms,dummy);

	//Local best fitness: lbfit[i][j] = i-th swarm, j-th individual
	vector< vector <double> > lbfit(nswarms,dummy2);
	//Local best chromosome: lbX[i][j][k] =  i-th swarm, j-th individual, k-th phenotype
	vector< vector< vector<double> > > lbX(nswarms,dummy3);

	double vwidth;										//Width of the search space
    vector<double> MINV(D),MAXV(D);						//Maximum and minumum velocity allowed

	// Initialise the particles (class Individual) positions, their velocities and their fitness to that of the deme
	for ( int i = 0; i<nswarms; i++){
		for ( int j = 0; j<NP; j++ ){
				X[i][j]	=	deme[NP*i+j].getDecisionVector();
				V[i][j]	=	deme[NP*i+j].getVelocity();
				fit[i][j]	=	deme[NP*i+j].getFitness();
		}
	}
   // Initialise the minimum and maximum velocity
   for ( int i = 0; i<D; i++ ) {
		vwidth = (UB[i]-LB[i])/vcoeff;
		MINV[i] = -1.0*vwidth;
		MAXV[i] = vwidth;
	}

   for (int i=0;i<nswarms;i++){
	// Initialise the global and local bests
	gbX[i]=X[i][0];
	gbfit[i]=fit[i][0];

	lbX[i]=X[i];			//at the first generation the local best position is the particle position
	lbfit[i]=fit[i];		//same for the fitness

	for (int j = 1; j<NP; j++){		//the int j = 1 jumps the first member as it is already set as the best
		if (fit[i][j] < gbfit[i]){
			gbfit[i] = fit[i][j];
			gbX[i] = X[i][j];
		}
	}
   }

   // Main PSO loop
   for (int iter = 0; iter<generations; iter++){
	//loop through the swarms
	for (int i = 0;i<nswarms;i++){
		//1 - move the particles and check that velocity and positions are in allowed range
		for (int j = 0; j< NP; j++){
			for (int k = 0; k< D; k++){

				//new velocity
				V[i][j][k] = omega * V[i][j][k] + eta1 * Pk::nextDouble(rng) * (lbX[i][j][k] - X[i][j][k]) + eta2 * Pk::nextDouble(rng) * (gbX[i][k] - X[i][j][k]);

				//check that it is within the allowed velocity range
				if ( V[i][j][k] > MAXV[k] )
					V[i][j][k] = MAXV[k];

				else if ( V[i][j][k] < MINV[k] )
					V[i][j][k] = MINV[k];

				//update position
				X[i][j][k] = X[i][j][k] + V[i][j][k];

				if (X[i][j][k] < LB[k])
					X[i][j][k] = Pk::nextDouble(rng) * (UB[k] - LB[k]) + LB[k];
				else if (X[i][j][k] > UB[k])
					X[i][j][k] = Pk::nextDouble(rng) * (UB[k] - LB[k]) + LB[k];
			}

			//We evaluate the new individual fitness now as to be able to update immediately the global best
			//in case a better solution is found
			fit[i][j] = problem.objfun(X[i][j]);
			//update local and global best
			if (fit[i][j] < lbfit[i][j]){
				lbfit[i][j] = fit[i][j];	//local best
				lbX[i][j] = X[i][j];
				if(fit[i][j] < gbfit[i]){
					gbfit[i] = fit[i][j];	//global best
					gbX[i]	= X[i][j];
				}
			}
		} //End of loop on the population members
	} //End of loop on the swarms

	//exchanges two random elements from randomly selected swarms
	if (iter % (int)5 == 0){
		int sw1 = (int)(Pk::nextDouble(rng)*nswarms);		//select 1st swarm
		int sw2 = (int)(Pk::nextDouble(rng)*nswarms);
		do{										        //endless loop if nswarms<2
		  sw2 = (int)(Pk::nextDouble(rng)*nswarms);		//selects 2nd swarm different from the first
		} while (sw2 == sw1);

		int in1 = (int)(Pk::nextDouble(rng)*NP);	        //select 1st individual
		int in2;
		do{										        //endless loop if nswarms<2
		  in2 = (int)(Pk::nextDouble(rng)*NP);		    //selects 2nd individual
		} while (in2 == in1);
		//swap position
		dummy = X[sw1][in1];
		X[sw1][in1]=X[sw2][in2];
		X[sw2][in2] = dummy;

		//swap velocity
		dummy = V[sw1][in1];
		V[sw1][in1]=V[sw2][in2];
		V[sw2][in2] = dummy;

		//swap local bests lbX
		dummy = lbX[sw1][in1];
		lbX[sw1][in1]=lbX[sw2][in2];
		lbX[sw2][in2] = dummy;

		//swap local best lbfit;
		dummy[0] = lbfit[sw1][in1];
		lbfit[sw1][in1]=lbfit[sw2][in2];
		lbfit[sw2][in2] = dummy[0];
	}//end of swap
   } // end of main PSO loop

   //we end by constructing the object Population containing the final results
   Population popout;
   Individual temp;
   for (int i=0; i<nswarms;i++){
   for (int j=0; j<NP; j++){
	temp.setDecisionVector(lbX[i][j]);
	temp.setFitness(lbfit[i][j]);
	temp.setVelocity(V[i][j]);
	popout.addIndividual(temp);
   }
   }
   return popout;
}

