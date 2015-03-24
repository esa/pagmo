/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

// Test code for the robust meta-problem

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include "../src/pagmo.h"

namespace pagmo{ namespace problem {

/// A dummy problem that outputs its input (i.e. x == objfun(x))
class white_box: public base
{
	public:
		white_box(size_type x_dim);
		base_ptr clone() const;

	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
};

white_box::white_box(size_type x_dim): base((int)x_dim, 0, x_dim, 0, 0, 0.0)
{
}

void white_box::objfun_impl(fitness_vector & f, const decision_vector & x) const
{
	f = x;
}

base_ptr white_box::clone() const
{
	return base_ptr(new white_box(*this));
}

}}

using namespace pagmo;

rng_uint32 urng;

const double EPS = 1e-8;

// Test strategy:
// Use white_box problem to probe and check if the perturbation is as expected
int test_robust(unsigned int dim, unsigned int n_trials, double rho, unsigned int seed = 0)
{
	std::cout << "[START] Testing robust meta-problem using white_box (rho = " << rho << ")" << std::endl;

	problem::white_box prob_white_box(dim);
	problem::robust robust_white_box(prob_white_box, n_trials, rho, seed);

	const unsigned int num_points = 500;

	population points(prob_white_box, num_points, seed);
	
	for(unsigned int i = 0; i < num_points; i++){

		// There has to be some perturbation
		double total_perturbation = 0.0;
		const double acceptable_min = EPS;

		// The perturbation has to be a mix of positive and negative ones
		bool found_positive_delta = false;
		bool found_negative_delta = false;

		std::cout << "\tPoint #" << i << ": ";

		const decision_vector& input_x = points.get_individual(i).cur_x;
		fitness_vector output_x = robust_white_box.objfun(input_x);

		// White box should behave like a white box
		if(output_x.size() != input_x.size()){
			std::cout << "FAILED: Something is wrong with the white box!" << std::endl;
		}

		// All of the perturbed points must be in the rho area
		for(unsigned int j = 0; j < output_x.size(); j++){
			double perturbation = output_x[j] - input_x[j];
			if(fabs(perturbation) > rho){
				std::cout << "FAILED @Dimension [" << j << "]: Exceeded rho box!" << std::endl;
				return 1;
			}
			total_perturbation += fabs(perturbation);
			found_positive_delta |= (perturbation > EPS);
			found_negative_delta |= (perturbation < -EPS);
		}

		if(total_perturbation < acceptable_min){
			std::cout << "FAILED: Perturbation is suspiciously small! (" << total_perturbation << ")" << std::endl;
			return 1;
		}

		if(!found_positive_delta){
			std::cout << "FAILED: No positive perturbation at all, too strange!" << std::endl;
			return 1;
		}

		if(!found_negative_delta){
			std::cout << "FAILED: No negative perturbation at all, too strange!" << std::endl;
			return 1;
		}

		// std::cout << "total_perturbation: " << total_perturbation << "   ";
		
		std::cout << "PASSED."<< std::endl;
	}

	std::cout << "[PASSED] Testing robust meta-problem using white_box (rho = " << rho << ")" << std::endl;

	return 0;
}

int main()
{
	return test_robust(10, 1, 0.001) ||
		   test_robust(20, 1, 0.01) ||
		   test_robust(30, 1, 0.1) ||
		   test_robust(40, 5, 0.5);
}
