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

// 16/05/2008: Initial version by Dario Izzo.

#include <iostream>
#include <vector>

#include "../../exceptions.h"
#include "../problems/GOproblem.h"
#include "PSO.h"

PSOalgorithm::PSOalgorithm(int generationsInit, const double &omegaInit, const double &eta1Init, const double &eta2Init,
	const double &vcoeffInit):generations(generationsInit),omega(omegaInit),eta1(eta1Init),eta2(eta2Init),vcoeff(vcoeffInit)
{
	if (generationsInit < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

Population PSOalgorithm::evolve(const Population &deme) const
{
    const GOProblem &problem = deme.problem();
    const std::vector<double>& LB = problem.getLB();
    const std::vector<double>& UB = problem.getUB();

	int n = deme.size();
	int m = LB.size();

	std::vector<double> dummy(m,0);							//used for initialisation purposes
	std::vector<std::vector<double> > X(n,dummy);
	std::vector<std::vector<double> > V(n,dummy);

	std::vector<double> fit(n);							//particle fitness

	double gbfit;										//global best fitness
	std::vector<double> gbX(m);								//global best chromosome

	std::vector<double> lbfit(n);							//local best fitness
	std::vector<std::vector<double> > lbX(n,dummy);				//local best chromosome

	double vwidth;										//Width of the search space
    std::vector<double> MINV(m),MAXV(m);						//Maximum and minumum velocity allowed

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
   for (size_t j = 0; j < generations; ++j){


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
   Population popout(problem,0);
   for (int i=0; i<n; i++){
	popout.push_back(Individual(problem,lbX[i],V[i]));
   }
   return popout;
}

void PSOalgorithm::log(std::ostream &s) const
{
	s << "PSO - generations:" << generations << " omega:" << omega << " eta1:" << eta1
		<< " eta2:" << eta2 << " vcoeff:" << vcoeff;
}
