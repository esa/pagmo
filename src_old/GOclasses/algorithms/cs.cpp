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

#include "../basic/population.h"
#include "../problems/base.h"
#include "base.h"
#include "cs.h"

namespace pagmo
{
namespace algorithm {

cs::cs(const double &minRange_):base(),range(0.25),
		reduxCoeff(0.5),minRange(minRange_)

{
	if (minRange_ > 1 || minRange_ <= 0 || minRange_>0.25) {
		pagmo_throw(value_error,"the minimum range must be smaller than one, positive and smaller than the starting range, (o44portebat studuisse)!!");
	}
}

cs::cs(const double &range_, const double &reduxCoeff_, const double &minRange_):base(),range(range_),
		reduxCoeff(reduxCoeff_),minRange(minRange_)

{
	if (reduxCoeff_ >= 1 || reduxCoeff_ <=0) {
		pagmo_throw(value_error,"the reduction coefficient must be smaller than one and positive, You Fool!!");
	}
	if (range_ > 1 || range_ <= 0) {
		pagmo_throw(value_error,"the starting range must be smaller than one and positive, You Fool!!");
	}
	if (minRange_ > 1 || minRange_ <= 0 || minRange_>range_) {
		pagmo_throw(value_error,"the minimum range must be smaller than one, positive and smaller than the starting range, (o44portebat studuisse)!!");
	}
}

population cs::evolve(const population &popin) const
{
	const problem::base &problem = popin.problem();
	if (popin.size() == 0) {
		return population(problem,0);
	}

	population retval(popin);
	const individual &startindividual = retval.extractBestIndividual();
	std::vector <double> x,newx;
	x = startindividual.get_decision_vector();
	double f = startindividual.get_fitness();
	double newf;
	size_t D = problem.getDimension();
	bool flag = false;
	std::vector <double> UB = problem.get_ub();
	std::vector <double> LB = problem.get_lb();

	double newrange=range;

	while (newrange > minRange) {
		flag = false;
		for (unsigned int i=0; i<D; i++) {
			newx=x;


			newx[i] = x[i] + newrange * (UB[i]-LB[i]);
			//feasibility correction
			if (newx[i] > UB [i]) newx[i]=UB[i];

			newf = problem.objfun(newx);
			if (newf < f) {
				f = newf;
				x = newx;
				flag=true;
				break; //accept
			}

			newx[i] = x[i] - newrange * (UB[i]-LB[i]);
			//feasibility correction
			if (newx[i] < LB [i]) newx[i]=LB[i];

			newf = problem.objfun(newx);
			if (newf < f) {  //accept
				f = newf;
				x = newx;
				flag=true;
				break;
			}
		}
		if (!flag) {
			newrange *= reduxCoeff;
		}
	} //end while

	retval.replace_best(individual(x, startindividual.get_velocity(), f));
	return retval;
}

void cs::log(std::ostream &s) const
{
	s << "CS - Starting Range:" << range << " FInal Range:" << minRange << " Reduction Coefficient:" << reduxCoeff;
}

}
}
