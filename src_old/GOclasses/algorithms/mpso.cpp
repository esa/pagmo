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

// 23/10/2008: Initial version by Dario Izzo.

#include <iostream>
#include <vector>

#include "../../exceptions.h"
#include "../basic/individual.h"
#include "../basic/population.h"
#include "../problems/base.h"
#include "base.h"
#include "mpso.h"

namespace pagmo
{
namespace algorithm {

mpso::mpso(int generationsInit, const double &omegaInit, const double &eta1Init, const double &eta2Init,
           const double &vcoeffInit, int nswarmsInit):base(),generations(generationsInit),omega(omegaInit),eta1(eta1Init),eta2(eta2Init),
		vcoeff(vcoeffInit),nswarms(nswarmsInit)
{
	if (generationsInit <= 0) {
		pagmo_throw(value_error,"number of generations must be positive");
	}
	if (nswarmsInit <= 0) {
		pagmo_throw(value_error,"number of swarms must be positive");
	}
}

population mpso::evolve(const population &deme) const
{
	const problem::base &problem = deme.problem();
	const std::vector<double>& LB = problem.get_lb();
	const std::vector<double>& UB = problem.get_ub();

	if (deme.size() % nswarms != 0) {
		pagmo_throw(value_error,"population size must be a multiple of the number of swarms");
	}

	const size_t NP = deme.size()/nswarms;		//potentially dangerous when the deme size is not divisible by the numberof swamrs
	const size_t D = LB.size();

	std::vector<double> dummy(D,0);								//used for initialisation purposes
	std::vector<double> dummy2(NP,0);							//used for initialisation purposes
	std::vector<std::vector<double> > dummy3(NP,dummy);				//used for initialisation purposes

	//Particle position: X[i][j][k] = i-th swarm, j-th individual, k-th phenotype
	std::vector<std::vector<std::vector<double> > > X(nswarms,dummy3);
	//Particle velocity: V[i][j][k] = i-th swarm, j-th individual, k-th phenotype
	std::vector<std::vector<std::vector<double> > > V(nswarms,dummy3);
	//Particle fitness: fit[i][j] = i-th swarm, j-th individual
	std::vector<std::vector<double> > fit(nswarms,dummy2);
	//Global best swarm fitness: gbfit[i] = i-th swarm
	std::vector<double> gbfit(nswarms,0);
	//Global best swarm individual: gbX[i][j] = i-th swarm, j-th phenotype
	std::vector<std::vector<double> > gbX(nswarms,dummy);

	//Local best fitness: lbfit[i][j] = i-th swarm, j-th individual
	std::vector<std::vector<double> > lbfit(nswarms,dummy2);
	//Local best chromosome: lbX[i][j][k] =  i-th swarm, j-th individual, k-th phenotype
	std::vector<std::vector<std::vector<double> > > lbX(nswarms,dummy3);

	double vwidth;										//Width of the search space
	std::vector<double> MINV(D),MAXV(D);						//Maximum and minumum velocity allowed

	// Initialise the particles (class individual) positions, their velocities and their fitness to that of the deme
	for (size_t i = 0; i < nswarms; ++i) {
		for (size_t j = 0; j < NP; ++j) {
			X[i][j]	= deme[NP*i+j].get_decision_vector();
			V[i][j]	= deme[NP*i+j].get_velocity();
			fit[i][j] = deme[NP*i+j].get_fitness();
		}
	}
	// Initialise the minimum and maximum velocity
	for (size_t i = 0; i < D; ++i) {
		vwidth = (UB[i]-LB[i])/vcoeff;
		MINV[i] = -1.0*vwidth;
		MAXV[i] = vwidth;
	}

	for (size_t i = 0; i < nswarms; ++i) {
		// Initialise the global and local bests
		gbX[i]=X[i][0];
		gbfit[i]=fit[i][0];

		lbX[i]=X[i];			//at the first generation the local best position is the particle position
		lbfit[i]=fit[i];		//same for the fitness

		for (size_t j = 1; j < NP; ++j) {		//the int j = 1 jumps the first member as it is already set as the best
			if (fit[i][j] < gbfit[i]) {
				gbfit[i] = fit[i][j];
				gbX[i] = X[i][j];
			}
		}
	}

	// Main PSO loop
	for (size_t iter = 0; iter < generations; ++iter) {
		//loop through the swarms
		for (size_t i = 0; i < nswarms; ++i) {
			//1 - move the particles and check that velocity and positions are in allowed range
			for (size_t j = 0; j < NP; ++j) {
				for (size_t k = 0; k< D; ++k) {

					//new velocity
					V[i][j][k] = omega * V[i][j][k] + eta1 * drng() * (lbX[i][j][k] - X[i][j][k]) + eta2 * drng() * (gbX[i][k] - X[i][j][k]);

					//check that it is within the allowed velocity range
					if ( V[i][j][k] > MAXV[k] )
						V[i][j][k] = MAXV[k];

					else if ( V[i][j][k] < MINV[k] )
						V[i][j][k] = MINV[k];

					//update position
					X[i][j][k] = X[i][j][k] + V[i][j][k];

					if (X[i][j][k] < LB[k])
						X[i][j][k] = drng() * (UB[k] - LB[k]) + LB[k];
					else if (X[i][j][k] > UB[k])
						X[i][j][k] = drng() * (UB[k] - LB[k]) + LB[k];
				}

				//We evaluate the new individual fitness now as to be able to update immediately the global best
				//in case a better solution is found
				fit[i][j] = problem.objfun(X[i][j]);
				//update local and global best
				if (fit[i][j] < lbfit[i][j]) {
					lbfit[i][j] = fit[i][j];	//local best
					lbX[i][j] = X[i][j];
					if (fit[i][j] < gbfit[i]) {
						gbfit[i] = fit[i][j];	//global best
						gbX[i]	= X[i][j];
					}
				}
			} //End of loop on the population members
		} //End of loop on the swarms

		//exchanges two random elements from randomly selected swarms
		if (iter % (int)5 == 0) {
			int sw1 = (int)(drng()*nswarms);		//select 1st swarm
			int sw2 = (int)(drng()*nswarms);
			do {										       //endless loop if nswarms<2
				sw2 = (int)(drng()*nswarms);		//selects 2nd swarm different from the first
			} while (sw2 == sw1);

			int in1 = (int)(drng()*NP);	        //select 1st individual
			int in2;
			do {										       //endless loop if nswarms<2
				in2 = (int)(drng()*NP);		    //selects 2nd individual
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

	//we end by constructing the object population containing the final results
	population popout(problem,0);
	for (size_t i = 0; i < nswarms; ++i) {
		for (size_t j = 0; j < NP; ++j) {
			popout.push_back(individual(lbX[i][j],V[i][j],lbfit[i][j]));
		}
	}
	return popout;
}

void mpso::log(std::ostream &s) const
{
	s << "MPSO - generations:" << generations << " omega:" << omega << " eta1:" << eta1
	<< " eta2:" << eta2 << " vcoeff:" << vcoeff
	<< " nswarms:" << nswarms;
}

}
}
