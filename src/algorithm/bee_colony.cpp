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
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>


namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations.
 * @param[in] onlooker_fraction fraction of onlooker bees
 * @param[in] scouts_number number of scout bees 
 * @param[in] limit number of tries a source of food is dropped if no better mutant solution is found
 * @param[in] phi define the random range where looking for neighbours  
 * @throws value_error if onlooker fraction is not in the [0,1] interval and ///scout_number > population_size
 */
bee_colony::bee_colony(int gen, double onlooker_fraction, int scouts_number, int limit, double phi):base(),m_gen(gen),m_onlooker_fraction(onlooker_fraction),m_scouts_number(scouts_number), m_limit(limit), m_phi(phi) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}

	if (onlooker_fraction < 0 || onlooker_fraction > 1) {
		pagmo_throw(value_error,"the fraction of onlooker bees must be in the [0,1] range");
	}

}

/// Clone method.
base_ptr bee_colony::clone() const
{
	return base_ptr(new bee_colony(*this));
}

double fitness_magnitude(fitness_vector fit);

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
	const population::size_type n_employers = (int) (pop.size() * (1-m_onlooker_fraction)); //number of employed bees
	const population::size_type n_onlookers = pop.size() - n_employers; //number of onlooker bees
	
	population::size_type scouts_number = m_scouts_number;
	//number of scout bees must be in the range [1,n_employers]
	if (scouts_number <=0) scouts_number = 1;
	if (scouts_number >= n_employers) scouts_number = n_employers;

	// Get out if there is nothing to do.
	if (n_employers == 0 || m_gen == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	std::vector<double> dummy(D,0);				//used for initialisation purposes
	std::vector<decision_vector > X(n_employers,dummy);	//set of food sources
	std::vector<fitness_vector> fit(n_employers);		//food sources fitness

	decision_vector temp_solution(D,0);

	std::vector<int> trial(n_employers,0);


	// Copy the particle positions, their velocities and their fitness
	for ( population::size_type i = 0; i<n_employers; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;		
	}

	double r = 0;
	// Main ABC loop
	for (int j = 0; j < m_gen; ++j) {
		//Send employed bees
		for (population::size_type ii = 0; ii< n_employers; ii++) {
			r = m_drng();
			int param2change = r*D;  //randomly determine the parameter to change

			r = m_drng();
			population::size_type neighbour = (int) (r*n_employers); //randomly chose a solution to be used to produce a mutant solution of solution ii

			while(neighbour == ii) { //randomly selected solution must be different from ii
				r = m_drng();
				neighbour = (int) (r*n_employers);
			}

			
			for(population::size_type i=0; i<D; i++) {
				temp_solution[i] = X[ii][i]; //copy local solution into temp_solution
			}

			r = m_drng();
			double selected_phi = m_phi*(r-0.5)*2; //randomly chose a phi value

			temp_solution[param2change] = X[ii][param2change] + selected_phi * (X[ii][param2change] - X[neighbour][param2change]);

			/*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        		if ((temp_solution[param2change])<lb[param2change]) 
           			temp_solution[param2change] = lb[param2change];
        		if ((temp_solution[param2change])>ub[param2change]) 
          			temp_solution[param2change] = ub[param2change];	

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
     		maxfit = fitness_magnitude(fit[0]);
  		for (population::size_type jj=1; jj<n_employers; jj++)
        	{
           		if (fitness_magnitude(fit[jj])>maxfit)
           		maxfit = fitness_magnitude(fit[jj]);
        	}

		//Calculate probability for an onlooker bee to chose a food source 
		double probability[n_employers]; 

		for (population::size_type jj=0; jj<n_employers; jj++)
        	{
         		probability[jj] = ( 0.9* (fitness_magnitude(fit[jj]) / maxfit) ) +0.1;
        	}
		
		population::size_type t = 0;
		population::size_type ii= 0;
		while(t < n_onlookers) {
			r = m_drng();
			if (r<probability[ii]) { //chose a food source depending on its probability to bee chosen
				t++;
				r = m_drng();
				int param2change = r*D;  //randomly determine the parameter to change

				r = m_drng();
				population::size_type neighbour = (int) (r*n_employers); //randomly chose a solution to be used to produce a mutant solution of solution ii

				while(neighbour == ii) { //randomly selected solution must be different from ii
					r = m_drng();
					neighbour = (int) (r*n_employers);
				}
				
				for(population::size_type i=0; i<D; i++) {
					temp_solution[i] = X[ii][i]; //copy local solution into temp_solution
				}

				r = m_drng();
				double selected_phi = m_phi*(r-0.5)*2; //randomly chose a phi value

				temp_solution[param2change] = X[ii][param2change] + selected_phi * (X[ii][param2change] - X[neighbour][param2change]);

				/*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
				if ((temp_solution[param2change])<lb[param2change]) 
					temp_solution[param2change] = lb[param2change];
				if ((temp_solution[param2change])>ub[param2change]) 
					temp_solution[param2change] = ub[param2change];	

				//if the new solution is better than the old one replace it with the mutant one and reset its trial counter 
				if(prob.compare_fitness(prob.objfun(temp_solution), fit[ii])) {
					X[ii][param2change] = temp_solution[param2change];
					pop.set_x(ii,X[ii]);
					prob.objfun(fit[ii], X[ii]); //update the fitness vector
				trial[ii] = 0;
				}
				else {
					trial[ii]++; //if the solution can't be improved incrase its the trial counter
				}
			}
			ii++;
			if (ii==n_employers-1) ii=0;
		}

		//Send scout bees
		for(population::size_type t=0; t < scouts_number; t++) {
			int maxtrialindex = 0;
			for (population::size_type ii=1; ii<n_employers; ii++)
			{
				 if (trial[ii] > trial[maxtrialindex])
				 	maxtrialindex = ii;
			}
			if(trial[maxtrialindex] >= m_limit)
			{
				//select a new random solution
				for(problem::base::size_type jj = 0; jj < D; jj++) {
					r = m_drng();
					X[maxtrialindex][jj] = r*(ub[jj]-lb[jj]) + lb[jj];
				}
				trial[maxtrialindex] = 0;
				pop.set_x(maxtrialindex,X[maxtrialindex]);
			}
		}




	} // end of main ABC loop
	
}

double fitness_magnitude(fitness_vector fit) {
	int size = sizeof(fit)/sizeof(double);
	double magnitude = 0;
	for(int i=0; i < size; i++) {
		magnitude += fit[i]*fit[i];
	}
	return sqrt(magnitude);
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
	s << "gen:" << m_gen << ' ';
	s << "onlooker_fraction:" << m_onlooker_fraction << ' ';
	s << "scouts_number:" << m_scouts_number << ' ';
	s << "limit:" << m_limit << ' ';
	s << "phi:" << m_phi << ' ';
	return s.str();
}

}} //namespaces
