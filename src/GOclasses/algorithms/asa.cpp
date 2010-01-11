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

#include <cmath>
#include <iostream>

#include "../../exceptions.h"
#include "../basic/individual.h"
#include "../basic/population.h"
#include "../problems/base.h"
#include "asa.h"

namespace pagmo
{
namespace algorithm {

asa::asa(int niterTotInit, const double &Ts_, const double &Tf_):base(),niterTot(niterTotInit),
		niterTemp(1),niterRange(20),Ts(Ts_),Tf(Tf_),StartStep(1)
{
	if (niterTotInit < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (Ts_ <= 0 || Tf_ <= 0 || Ts_ <= Tf_) {
		pagmo_throw(value_error,"temperatures must be positive and Ts must be greater than Tf");
	}
}

asa::asa(int niterTotInit, const double &Ts_, const double &Tf_, const int niterTemp_, const int niterRange_, const double StartStep_):
		base(),niterTot(niterTotInit),niterTemp(niterTemp_),niterRange(niterRange_),Ts(Ts_),Tf(Tf_),StartStep(StartStep_)
{
	if (niterTotInit < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (Ts_ <= 0 || Tf_ <= 0 || Ts_ <= Tf_) {
		pagmo_throw(value_error,"temperatures must be positive and Ts must be greater than Tf");
	}
	if (niterTemp_ < 0) {
		pagmo_throw(value_error,"number of iteration before adjusting the temperature must be positive");
	}
	if (niterRange_ < 0) {
		pagmo_throw(value_error,"number of iteration before adjusting the neighbourhood must be positive");
	}
	if (StartStep_ < 0 || StartStep_>1) {
		pagmo_throw(value_error,"Initial range must be between 0 and 1");
	}
}

population asa::evolve(const population &pop) const
{
	const problem::base &problem = pop.problem();
	if (pop.size() == 0) {
		return population(problem,0);
	}
	const std::vector<double> &LB = problem.get_lb();
	const std::vector<double> &UB = problem.get_ub();
	const size_t SolDim = LB.size();
	if (SolDim == 0) {
		pagmo_throw(value_error,"problem's size cannot be null");
	}
	const size_t niterOuter = niterTot / (niterTemp * niterRange * SolDim);
	if (niterOuter == 0) {
		pagmo_throw(value_error,"niterOuter cannot be null");
	}
	population retval(pop);
	const individual &x0 = retval.extractBestIndividual();
	const double Tcoeff = std::pow(Tf/Ts,1.0/(double)(niterOuter));
	std::vector<double> xNEW = x0.get_decision_vector(), xOLD = xNEW;
	double fNEW = x0.get_fitness(), fOLD = fNEW;
	if (xNEW.size() != SolDim) {
		pagmo_throw(value_error,"discrepancy between individual size and problem size.");
	}
	std::vector<double> Step(SolDim,StartStep);
	std::vector<int> acp(SolDim,0) ;
	double ratio = 0, currentT = Ts, prob = 0,  r = 0;

	//Main SA loops
	for (size_t jter = 0; jter < niterOuter; ++jter) {
		for (size_t mter = 0; mter < niterTemp; ++mter) {
			for (size_t kter = 0; kter < niterRange; ++kter) {
				size_t nter =  (size_t)(drng() * SolDim);
				for (size_t numb = 0; numb < SolDim ; ++numb) {
					nter = (nter + 1) % SolDim;
					//We modify the current point actsol by mutating its nter component within
					//a Step that we will later adapt
					r = 2.0 * drng() - 1.0; //random number in [-1,1]
					xNEW[nter] = xOLD[nter] + r * Step[nter] * ( UB[nter] - LB[nter] );

					// If new solution produced is infeasible ignore it
					if ((xNEW[nter] > UB[nter]) || (xNEW[nter] < LB[nter])) {
						xNEW[nter]=xOLD[nter];
						continue;
					}
					//And we valuate the objective function for the new point
					fNEW = problem.objfun(xNEW);

					// We decide wether to accept or discard the point
					if (fNEW < fOLD) {
						//accept
						xOLD[nter] = xNEW[nter];
						fOLD = fNEW;
						acp[nter]++;	//Increase the number of accepted values
					} else {
						//test it with Boltzmann to decide the acceptance
						prob = exp ( (fOLD - fNEW )/ currentT );

						// we compare prob with a random probability.
						if (prob > drng()) {
							xOLD[nter] = xNEW[nter];
							fOLD = fNEW;
							acp[nter]++;	//Increase the number of accepted values
						} else {
							xNEW[nter] = xOLD[nter];
						}
					} // end if
				} // end for(nter = 0; ...
			} // end for(kter = 0; ...


			// adjust the step (adaptively)
			for (size_t iter = 0; iter < SolDim; ++iter) {
				ratio = (double)acp[iter]/(double)niterRange;
				acp[iter] = 0;  //reset the counter
				if (ratio > .6) {
					//too many acceptances, increase the step by a factor 3 maximum
					Step[iter] = Step [iter] * (1 + 2 *(ratio - .6)/.4);
				} else {
					if (ratio < .4) {
						//too few acceptance, decrease the step by a factor 3 maximum
						Step [iter]= Step [iter] / (1 + 2 * ((.4 - ratio)/.4));
					};
				};
				//And if it becomes too large, reset it to its initial value
				if ( Step[iter] > StartStep ) {
					Step [iter] = StartStep;
				};
			}
		}
		// Cooling schedule
		currentT *= Tcoeff;
	}

	if (fOLD < x0.get_fitness()) {
		retval.replace_best(individual(xOLD,x0.get_velocity(),fOLD));
	}
	return retval;
}

void asa::log(std::ostream &s) const
{
	s << "ASA - Iter:" << niterTot << " Ts:" << Ts << " Tf:" << Tf
	<< " niterTemp:" << niterTemp << " niterRange:" << niterRange
	<< " StartStep:" << StartStep;
}

}
} // Close namespaces.
