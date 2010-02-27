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

// 18/05/2008: Initial version by Dario Izzo.

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../../exceptions.h"
#include "../basic/individual.h"
#include "../basic/population.h"
#include "../problems/base.h"
#include "base.h"
#include "de.h"

namespace pagmo
{
namespace algorithm {

de::de(int g, const double &F_, const double &CR_, int s):base(),generations(g),F(F_),CR(CR_),strategy(s)
{
	if (g < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (s < 1 || s > 10) {
		pagmo_throw(value_error,"strategy index must be between 1 and 10");
	}
	if (CR <= 0 || F <= 0 || CR >= 1 || F >= 1) {
		pagmo_throw(value_error,"the F and CR parameters must be in the ]0,1[ range");
	}
}

population de::evolve(const population &deme) const
{
	const problem::base &problem = deme.problem();
	const std::vector<double> &LB = problem.get_lb();
	const std::vector<double> &UB = problem.get_ub();
	const size_t D = LB.size(), NP = deme.size();
	// TODO: this needs to go away.
	if (NP < 6) {
		pagmo_throw(value_error,"for DE at least 6 individuals in the population are needed");
	}

	std::vector<double> dummy(D), tmp(D);						//dummy is used for initialisation purposes, tmp to contain
	//the mutated candidate
	std::vector<std::vector<double> > popold(NP,dummy), popnew(NP,dummy), popswap(NP,dummy);
	std::vector<double> fit(NP);								//chromosome fitness

	double newfitness;									//new fitness of the mutaded candidate
	double gbfit;										//global best fitness
	std::vector<double> gbX(D);								//global best chromosome
	std::vector<double> gbIter(D);							//best chromosome of current iteration

	// Initialise the chromosome (class individual) values, and their fitness to that of the deme
	for (size_t i = 0; i < NP; ++i) {
		popold[i] = deme[i].get_decision_vector();
		popnew[i] = deme[i].get_decision_vector();
		fit[i] = deme[i].get_fitness();
	}

	// Initialise the global bests
	gbX=popold[0];
	gbfit=fit[0];

	for (size_t i = 1; i < NP; ++i) {		//the int i = 1 jumps the first member as it is already set as the best
		if (fit[i] < gbfit) {
			gbfit = fit[i];			// save best member ever
			gbX = popold[i];
		}
	}
	gbIter = gbX;				// save best member of generation

	// Main DE iterations
	size_t r1,r2,r3,r4,r5;	//indexes to the selected population members
	for (size_t gen = 0; gen < generations; ++gen) {
		//Start of the loop through the deme
		for (size_t i = 0; i < NP; ++i) {
			do {                       /* Pick a random population member */
				/* Endless loop for NP < 2 !!!     */
				r1 = (size_t)(drng()*NP);
			} while (r1==i);

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 3 !!!     */
				r2 = (size_t)(drng()*NP);
			} while ((r2==i) || (r2==r1));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 4 !!!     */
				r3 = (size_t)(drng()*NP);
			} while ((r3==i) || (r3==r1) || (r3==r2));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 5 !!!     */
				r4 = (size_t)(drng()*NP);
			} while ((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 6 !!!     */
				r5 = (size_t)(drng()*NP);
			} while ((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));


			/*-------DE/best/1/exp--------------------------------------------------------------------*/
			/*-------Our oldest strategy but still not bad. However, we have found several------------*/
			/*-------optimization problems where misconvergence occurs.-------------------------------*/
			if (strategy == 1) { /* strategy DE0 (not in our paper) */
				tmp = popold[i];
				size_t n = (size_t)(drng()*D), L = 0;
				do {
					tmp[n] = gbIter[n] + F*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%D;
					++L;
				} while ((drng() < CR) && (L < D));
			}

			/*-------DE/rand/1/exp-------------------------------------------------------------------*/
			/*-------This is one of my favourite strategies. It works especially well when the-------*/
			/*-------"gbIter[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5---------*/
			/*-------as a first guess.---------------------------------------------------------------*/
			else if (strategy == 2) { /* strategy DE1 in the techreport */
				tmp = popold[i];
				size_t n = (size_t)(drng()*D), L = 0;
				do {
					tmp[n] = popold[r1][n] + F*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%D;
					++L;
				} while ((drng() < CR) && (L < D));
			}

			/*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
			/*-------This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.------*/
			/*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
			/*-------should play around with all three control variables.----------------------------*/
			else if (strategy == 3) { /* similiar to DE2 but generally better */
				tmp = popold[i];
				size_t n = (size_t)(drng()*D), L = 0;
				do {
					tmp[n] = tmp[n] + F*(gbIter[n] - tmp[n]) + F*(popold[r1][n]-popold[r2][n]);
					n = (n+1)%D;
					++L;
				} while ((drng() < CR) && (L < D));
			}
			/*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
			else if (strategy == 4) {
				tmp = popold[i];
				size_t n = (size_t)(drng()*D), L = 0;
				do {
					tmp[n] = gbIter[n] +
					         (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					n = (n+1)%D;
					++L;
				} while ((drng() < CR) && (L < D));
			}
			/*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
			else if (strategy == 5) {
				tmp = popold[i];
				size_t n = (size_t)(drng()*D), L = 0;
				do {
					tmp[n] = popold[r5][n] +
					         (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					n = (n+1)%D;
					++L;
				} while ((drng() < CR) && (L < D));
			}

			/*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

			/*-------DE/best/1/bin--------------------------------------------------------------------*/
			else if (strategy == 6) {
				tmp = popold[i];
				size_t n = (size_t)(drng()*D);
				for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
					if ((drng() < CR) || L + 1 == D) { /* change at least one parameter */
						tmp[n] = gbIter[n] + F*(popold[r2][n]-popold[r3][n]);
					}
					n = (n+1)%D;
				}
			}
			/*-------DE/rand/1/bin-------------------------------------------------------------------*/
			else if (strategy == 7) {
				tmp = popold[i];
				size_t n = (size_t)(drng()*D);
				for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
					if ((drng() < CR) || L + 1 == D) { /* change at least one parameter */
						tmp[n] = popold[r1][n] + F*(popold[r2][n]-popold[r3][n]);
					}
					n = (n+1)%D;
				}
			}
			/*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
			else if (strategy == 8) {
				tmp = popold[i];
				size_t n = (size_t)(drng()*D);
				for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
					if ((drng() < CR) || L + 1 == D) { /* change at least one parameter */
						tmp[n] = tmp[n] + F*(gbIter[n] - tmp[n]) + F*(popold[r1][n]-popold[r2][n]);
					}
					n = (n+1)%D;
				}
			}
			/*-------DE/best/2/bin--------------------------------------------------------------------*/
			else if (strategy == 9) {
				tmp = popold[i];
				size_t n = (size_t)(drng()*D);
				for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
					if ((drng() < CR) || L + 1 == D) { /* change at least one parameter */
						tmp[n] = gbIter[n] +
						         (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					}
					n = (n+1)%D;
				}
			}
			/*-------DE/rand/2/bin--------------------------------------------------------------------*/
			else if (strategy == 10) {
				tmp = popold[i];
				size_t n = (size_t)(drng()*D);
				for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
					if ((drng() < CR) || L + 1 == D) { /* change at least one parameter */
						tmp[n] = popold[r5][n] +
						         (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					}
					n = (n+1)%D;
				}
			}


			/*=======Trial mutation now in tmp[]. Force feasibility and how good this choice really was.==================*/
			// a) feasibility
			size_t i2 = 0;
			while (i2<D) {
				if ((tmp[i2] < LB[i2]) || (tmp[i2] > UB[i2]))
					tmp[i2] = drng()*(UB[i2]-LB[i2]) + LB[i2];
				++i2;
			}
			newfitness = problem.objfun(tmp);    /* Evaluate new vector in tmp[] */
			//cout << "\n" << newfitness << endl;
			//for (int g=0;g<NP;g++){
			//cout << tmp[g] << " ";
			//}
			//b) how good?
			if (newfitness <= fit[i]) {  /* improved objective function value ? */
				fit[i]=newfitness;
				popnew[i] = tmp;
				if (newfitness<gbfit) {        /* Was this a new minimum for the deme? */
					/* if so...*/
					gbfit=newfitness;          /* reset gbfit to new low...*/
					gbX=tmp;
				}
			} else {
				popnew[i] = popold[i];
			}
			/* swap population arrays. New generation becomes old one */

		}//End of the loop through the deme

		/* Save best population member of current iteration */
		gbIter = gbX;

		/* swap population arrays. New generation becomes old one */
		for (size_t i = 0; i < NP; ++i) {
			popswap[i] = popold[i];
			popold[i] = popnew[i];
			popnew[i] = popswap[i];
		}

	}//end main DE iterations

	//we end by constructing the object population containing the final results
	population popout(problem,0);
	std::vector<double> Xini(D),Vfin(D);
	for (size_t i = 0; i < NP; ++i) {
		Xini = deme[i].get_decision_vector();
		for (size_t j = 0; j < D; ++j) {
			Vfin[j] = popold[i][j] - Xini[j];
		}
		//X[i] - deme[i].get_decision_vector());  DOES NOT WORK AS VECTOR CLASS DOES NOT ACCEPT MINUS AS OPERATOR

		popout.push_back(individual(popold[i], deme[i].get_velocity(), fit[i]));
	}
	return popout;
}

void de::log(std::ostream &s) const
{
	s << "DE - generations:" << generations << " F:" << F << " CR:" << CR
	<< " strategy:" << strategy;
}

std::string de::id_object() const
{
	std::stringstream tmp;
	tmp << id_name() << "_" << strategy;
	return tmp.str();
}

}
}
