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

#include "bee_colony.h"
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>


namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] iter number of iterations.
 * @param[in] limit number of tries a source of food is dropped if no better mutant solution is found
 * @throws value_error if number of iterations or limit are negative
 */
bee_colony::bee_colony(int iter, int limit):base(),m_iter(iter), m_limit(limit) {
	if (iter < 0) {
		pagmo_throw(value_error,"number of iterations must be nonnegative");
	}

	if (limit < 0) {
		pagmo_throw(value_error,"limit value must be nonnegative");
	}

}

/// Clone method.
base_ptr bee_colony::clone() const
{
	return base_ptr(new bee_colony(*this));
}

/// Evolve implementation.
/**
 * Run the ABC algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void bee_colony::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = (int) pop.size();

	//We perform some checks to determine wether the problem/population are suitable for ABC
	if ( D - prob.get_i_dimension() == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for ABC to optimise");
	}

	if ( prob.get_f_dimension() != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and ABC is not suitable to solve it");
	}


	
	// Get out if there is nothing to do.
	if (NP == 0 || m_iter == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy);	//set of food sources
	std::vector<fitness_vector> fit(NP);		//food sources fitness

	decision_vector temp_solution(D,0);

	std::vector<int> trial(NP,0);

	std::vector<double> probability(NP);

	population::size_type neighbour = 0;

	int param2change = 0;

	double r = 0;

	// Copy the particle positions, their velocities and their fitness
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;		
	}

	// Main ABC loop
	for (int j = 0; j < m_iter; ++j) {
		//Send employed bees
		for (population::size_type ii = 0; ii< NP; ++ii) {
			param2change = boost::uniform_int<int>(0,D-1)(m_urng); 
			do{
				neighbour = boost::uniform_int<int>(0,NP-1)(m_urng); //randomly chose a solution to be used to produce a mutant solution of solution ii
			}
			while(neighbour == ii); //randomly selected solution must be different from ii

			
			for(population::size_type i=0; i<D; ++i) {
				temp_solution[i] = X[ii][i]; //copy local solution into temp_solution
			}

			temp_solution[param2change] = X[ii][param2change] + boost::uniform_real<double>(-1,1)(m_drng) * (X[ii][param2change] - X[neighbour][param2change]);

			/*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        		if ((temp_solution[param2change])<lb[param2change]) {
           			temp_solution[param2change] = lb[param2change];
			}
        		if ((temp_solution[param2change])>ub[param2change]) {
          			temp_solution[param2change] = ub[param2change];
			}

			//if the new solution is better than the old one replace it with the mutant one and reset its trial counter 
			if(prob.compare_fitness(prob.objfun(temp_solution), fit[ii])) {
				X[ii][param2change] = temp_solution[param2change];
				pop.set_x(ii,X[ii]);
				prob.objfun(fit[ii], X[ii]); //update the fitness vector
				trial[ii] = 0;
			}
			else {
				trial[ii]++; //if the solution can't be improved incrase its  trial counter
			}
		} //End of loop on the population members

		//Send onlooker bees
     		double maxfit;
     		maxfit = fit[0][0];
  		for (population::size_type jj=1; jj<NP; ++jj)
        	{
           		if (fit[jj][0]>maxfit)
           		maxfit = fit[jj][0];
        	}

		//Calculate probability for an onlooker bee to chose a food source 
		for (population::size_type jj=0; jj<NP; ++jj)
        	{
         		probability[jj] = ( 0.9* (fit[jj][0] / maxfit) ) +0.1;
        	}
		
		population::size_type t = 0;
		population::size_type ii= 0;
		while(t < NP) {
			r = m_drng();
			if (r<probability[ii]) { //chose a food source depending on its probability to bee chosen
				t++;
				param2change = boost::uniform_int<int>(0,D-1)(m_urng); 
				do{
					neighbour = boost::uniform_int<int>(0,NP-1)(m_urng); //randomly chose a solution to be used to produce a mutant solution of solution ii
				}
				while(neighbour == ii); //randomly selected solution must be different from ii

				
				for(population::size_type i=0; i<D; ++i) {
					temp_solution[i] = X[ii][i]; //copy local solution into temp_solution
				}

				temp_solution[param2change] = X[ii][param2change] + boost::uniform_real<double>(-1,1)(m_drng) * (X[ii][param2change] - X[neighbour][param2change]);

				/*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
				if ((temp_solution[param2change])<lb[param2change]) {
					temp_solution[param2change] = lb[param2change];
				}
				if ((temp_solution[param2change])>ub[param2change]) {
					temp_solution[param2change] = ub[param2change];
				}

				//if the new solution is better than the old one replace it with the mutant one and reset its trial counter 
				if(prob.compare_fitness(prob.objfun(temp_solution), fit[ii])) {
					X[ii][param2change] = temp_solution[param2change];
					pop.set_x(ii,X[ii]);
					prob.objfun(fit[ii], X[ii]); //update the fitness vector
					trial[ii] = 0;
				}
				else {
					trial[ii]++; //if the solution can't be improved incrase its  trial counter
				}
			}
			ii++;
			if (ii==NP-1) ii=0;
		}

		//Send scout bees
		int maxtrialindex = 0;
		for (population::size_type ii=1; ii<NP; ++ii)
		{
			 if (trial[ii] > trial[maxtrialindex]) {
				maxtrialindex = ii;
			 }
		}
		if(trial[maxtrialindex] >= m_limit)
		{
			//select a new random solution
			for(problem::base::size_type jj = 0; jj < D; ++jj) {
				X[maxtrialindex][jj] = boost::uniform_real<double>(lb[jj],ub[jj])(m_drng);
			}
			trial[maxtrialindex] = 0;
			pop.set_x(maxtrialindex,X[maxtrialindex]);
		}




	} // end of main ABC loop
	
}

/// Algorithm name
std::string bee_colony::get_name() const
{
	return "Artificial Bee Colony optimization";
}


/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string bee_colony::human_readable_extra() const
{
	std::ostringstream s;
	s << "iter:" << m_iter << ' ';
	s << "limit:" << m_limit << ' ';
	return s.str();
}

}} //namespaces
