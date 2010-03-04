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

// 16/05/2008: Initial version by Dario Izzo.

#include <iostream>
#include <string>
#include <vector>

#include "../../exceptions.h"
#include "../problems/base.h"
#include "base.h"
#include "pso.h"

namespace pagmo
{
namespace algorithm {

pso::pso(int generations_, const double &omega_, const double &eta1_, const double &eta2_,
         const double &vcoeff_):base(),generations(generations_),omega(omega_),eta1(eta1_),eta2(eta2_),vcoeff(vcoeff_), strategy(3)
{
	if (generations_ < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (vcoeff_ < 0 || vcoeff_ > 1) {
		pagmo_throw(value_error,"Velocity coeficient must be between 0 and 1");
	}
}

pso::pso(int generations_):base(),generations(generations_),omega(0.65),eta1(2.0),eta2(2.0),vcoeff(1), strategy(3)
{
	if (generations_ < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

pso::pso(int generations_, const double &omega_, const double &eta1_, const double &eta2_,
         const double &vcoeff_, const int &strategy_):base(),generations(generations_),omega(omega_),eta1(eta1_),eta2(eta2_),vcoeff(vcoeff_), strategy(strategy_)
{
	if (generations_ < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (vcoeff_ < 0 || vcoeff_ > 1) {
		pagmo_throw(value_error,"Velocity coeficient must be between 0 and 1");
	}
	if (strategy_ < 1|| strategy_ > 4) {
		pagmo_throw(value_error,"Staretgy for PSO must be an integer between 1 and 4");
	}
}

pso::pso(int generations_, const int &strategy_):base(),generations(generations_),omega(0.65),eta1(2.0),eta2(2.0),vcoeff(1), strategy(strategy_)
{
	if (generations_ < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (strategy_ < 1|| strategy_ > 4) {
		pagmo_throw(value_error,"Staretgy for PSO must be an integer between 1 and 4");
	}
}

population pso::evolve(const population &deme) const
{
	const problem::base &problem = deme.problem();
	const std::vector<double>& LB = problem.get_lb();
	const std::vector<double>& UB = problem.get_ub();

	int n = deme.size();
	int m = LB.size();

	std::vector<double> dummy(m,0);				//used for initialisation purposes
	std::vector<std::vector<double> > X(n,dummy);
	std::vector<std::vector<double> > V(n,dummy);

	std::vector<double> fit(n);				//particle fitness
	double gbfit;						//global best fitness
	std::vector<double> gbX(m);				//global best chromosome
	std::vector<double> lbfit(n);				//local best fitness
	std::vector<std::vector<double> > lbX(n,dummy);		//local best chromosome

	double vwidth;						//Width of the search space
	std::vector<double> MINV(m),MAXV(m);			//Maximum and minumum velocity allowed

	// Initialise the particle (class individual) positions, their velocities and their fitness to that of the deme
	for ( int i = 0; i<n; i++ ) {
		X[i]	=	deme[i].get_decision_vector();
		V[i]	=	deme[i].get_velocity();
		fit[i]	=	deme[i].get_fitness();
	}

	// Initialise the minimum and maximum velocity
	for ( int i = 0; i<m; i++ ) {
		vwidth = (UB[i]-LB[i]) * vcoeff;
		MINV[i] = -1.0*vwidth;
		MAXV[i] = vwidth;
	}

	// Initialise the global and local bests
	gbX=X[0];
	gbfit=fit[0];

	lbX=X;			//at the first generation the local best position is the particle position
	lbfit=fit;		//same for the fitness

	for (int i = 1; i<n; i++) {		//the int i = 1 jumps the first member as it is already set as the best
		if (fit[i] < gbfit) {
			gbfit = fit[i];
			gbX = X[i];
		}
	}


	double r1,r2 = 0;
	// Main PSO loop
	for (size_t j = 0; j < generations; ++j) {
		//For each particle in the swarm
		for (int ii = 0; ii< n; ii++) {

			/*-------PSO canonical--------------------------------------------------------------------*/
			/*-------The classical PSO velocity update startegy---------------------------------------*/
			if (strategy==1) {
				r1 = drng();
				r2 = drng();
				for (int jj = 0; jj< m; jj++) {
					V[ii][jj] = omega * V[ii][jj] + eta1 * r1 * (lbX[ii][jj] - X[ii][jj]) + eta2 * r2 * (gbX[jj] - X[ii][jj]);
				}
			}

			/*-------PSO canonical with equal random weights of social and cognitive components-------*/
			/*-------In our experience few problems benefit a lot from having r1=r2. You may check----*/
			/*-------with Rastrigin-------------------------------------------------------------------*/
			if (strategy==2) {
				r1 = drng();
				for (int jj = 0; jj< m; jj++) {
					V[ii][jj] = omega * V[ii][jj] + eta1 * r1 * (lbX[ii][jj] - X[ii][jj]) + eta2 * r1 * (gbX[jj] - X[ii][jj]);
				}
			}

			/*-------Variant of PSO strategy 1 that has r1 and r2 randomly generated for each----------*/
			/*-------component. This is also the version that was implemented--------------------------*/
			/*-------in early versions of pagmo--------------------------------------------------------*/
			if (strategy==3) {
				for (int jj = 0; jj< m; jj++) {
					r1 = drng();
					r2 = drng();
					V[ii][jj] = omega * V[ii][jj] + eta1 * r1 * (lbX[ii][jj] - X[ii][jj]) + eta2 * r2 * (gbX[jj] - X[ii][jj]);
				}
			}

			/*-------Variant of PSO strategy 2 that has r1 randomly generated for each-----------------*/
			/*-------component.------------------------------------------------------------------------*/
			if (strategy==4) {
				for (int jj = 0; jj< m; jj++) {
					r1 = drng();
					V[ii][jj] = omega * V[ii][jj] + eta1 * r1 * (lbX[ii][jj] - X[ii][jj]) + eta2 * r1 * (gbX[jj] - X[ii][jj]);
				}
			}

			//We now check that the velocity does not exceed the maximum allowed per component
			//and we perform the position update and the feasibility correction
			for (int jj = 0; jj< m; jj++) {

				if ( V[ii][jj] > MAXV[jj] )
					V[ii][jj] = MAXV[jj];

				else if ( V[ii][jj] < MINV[jj] )
					V[ii][jj] = MINV[jj];

				//update position
				X[ii][jj] = X[ii][jj] + V[ii][jj];

				//feasibility correction
				if (X[ii][jj] < LB[jj])
					X[ii][jj] = drng() * (UB[jj] - LB[jj]) + LB[jj];

				else if (X[ii][jj] > UB[jj])
					X[ii][jj] = drng() * (UB[jj] - LB[jj]) + LB[jj];
			}

			//We evaluate here the new individual fitness as to be able to update the global best in real time
			fit[ii] = problem.objfun(X[ii]);
			//update local and global best
			if (fit[ii] < lbfit[ii]) {
				lbfit[ii] = fit[ii];	//local best
				lbX[ii] = X[ii];
				if (fit[ii] < gbfit) {
					gbfit = fit[ii];	//global best
					gbX = X[ii];
				}
			}
		} //End of loop on the population members
	} // end of main PSO loop

	//we end by constructing the object population containing the final results
	population popout(problem,0);
	for (int i=0; i<n; i++) {
		popout.push_back(individual(lbX[i],V[i],lbfit[i]));
	}
	return popout;
}

void pso::log(std::ostream &s) const
{
	s << "PSO - generations:" << generations << " omega:" << omega << " eta1:" << eta1
	<< " eta2:" << eta2 << " vcoeff:" << vcoeff << " strategy:" << strategy;
}

std::string pso::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << strategy;
	return tmp.str();
}

}
}
