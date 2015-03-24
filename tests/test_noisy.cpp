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

// Test code for the noisy meta-problem

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include "../src/pagmo.h"
#include "test.h"

using namespace pagmo;

const double EPS = 10e-9;
rng_uint32 rng_seed_provider;

// Run the batch test by constructing noisy meta-problems against probs given the parameters.
// Test that the expected value of the noise is close to the mean after enough trials.
// noise_mean, noise_stddev: Params of the noise distribution
// num_trials: How many samples to take for calculating the expected value of the noise a some point.
// tol: Tolerance for the closeness
int test_noisy(const std::vector<problem::base_ptr> & probs, double noise_mean, double noise_stddev, int num_trials, double tol)
{
	
	std::cout << "Start batch testing with noise_mean = " << noise_mean << " and " 
	<< " noise_stddev = " << noise_stddev << ", sampling " << num_trials <<" times."<<std::endl;

	for(unsigned int i = 0; i < probs.size(); i++){

		int dim = probs[i]->get_dimension();

		pagmo::problem::noisy prob_noisy(*(probs[i]),
										 1,
										 noise_mean,
										 noise_stddev,
										 pagmo::problem::noisy::NORMAL,
										 i*177 + 23);

		std::cout<< std::setw(40) << prob_noisy.get_name();

		// Check that bounds are set correctly in the meta-problem
		if(!is_eq_vector(prob_noisy.get_lb(), probs[i]->get_lb(), EPS) ||
		   !is_eq_vector(prob_noisy.get_ub(), probs[i]->get_ub(), EPS)){
			std::cout<<" bounds failed!"<<std::endl;
			return 1;
		}
		std::cout << " bounds passed. ";

		// Generate a test point (middle of the search space)
		decision_vector x(dim);
		for(unsigned int xi = 0; xi < x.size(); xi++){
			x[xi] = (probs[i]->get_lb()[xi] + probs[i]->get_ub()[xi]) / 2.0;
		}

		// Check that the sample mean of the noise is the same as specified by noise_mean
		fitness_vector f_noiseless = probs[i]->objfun(x);
		fitness_vector noise_sample_mean(f_noiseless.size(), 0.0);
		fitness_vector noise_sample_variance(f_noiseless.size(), 0.0);

		// Compute the sample mean and variance of the generated noise
		// Variance sanity will be check after mean's sanity, thus directly minus the noise_mean
		// (won't get to checking variance is noise_mean is not ok)
		for(int tr = 0; tr < num_trials; tr++){
			prob_noisy.set_seed(rng_seed_provider()); // Automatically resets the cache
			fitness_vector f_cur = prob_noisy.objfun(x);
			for(unsigned int fi = 0; fi < f_cur.size(); fi++){
				noise_sample_mean[fi] += (f_cur[fi] - f_noiseless[fi]) / (double)num_trials;
				noise_sample_variance[fi] += (f_cur[fi] - f_noiseless[fi] - noise_mean) * (f_cur[fi] - f_noiseless[fi] - noise_mean) / (double)num_trials;
			}
		}

		// Check if mean is as specified
		for(unsigned int fi = 0; fi < noise_sample_mean.size(); fi++){
			if(fabs(noise_sample_mean[fi] - noise_mean) > tol){
				std::cout<< " mean of noise is strange!"<<std::endl;
				std::cout<< noise_sample_mean << ", expected = " << noise_mean << std::endl;
				return 1;
			}
		}
		std::cout << " mean of noise passed. ";

		// Check if variance is as specified
		for(unsigned int fi = 0; fi < noise_sample_variance.size(); fi++){
			if(fabs(noise_sample_variance[fi] - noise_stddev*noise_stddev) > tol){
				std::cout<< " variance of noise is strange!"<<std::endl;
				std::cout<< noise_sample_variance << ", expected = " << " " << noise_stddev*noise_stddev << std::endl;
				return 1;
			}
		}
		std::cout << " variance of noise passed. " << std::endl;

		// Can also check on constraints...
		// But currently both generated from the same RNG and distribution.
	}

	// All good
	return 0;
}


// Run the batch test by constructing noisy meta-problems against probs given the parameters.
// Test similar to above, except that now the two noise params are interpreted as the lower
// and upper bounds of the uniform noise. The sanity check is thus a bit different.
int test_noisy_uniform(const std::vector<problem::base_ptr> & probs, double noise_lb, double noise_ub, int num_trials, double tol)
{
	
	std::cout << "Start batch testing with noise_lb = " << noise_lb << " and " 
	<< " noise_ub = " << noise_ub << ", sampling " << num_trials <<" times."<<std::endl;

	for(unsigned int i = 0; i < probs.size(); i++){

		int dim = probs[i]->get_dimension();

		pagmo::problem::noisy prob_noisy(*(probs[i]),
										 1,
										 noise_lb,
										 noise_ub,
										 pagmo::problem::noisy::UNIFORM,
										 i*177 + 23);

		std::cout<< std::setw(40) << prob_noisy.get_name();

		// Check that bounds are set correctly in the meta-problem
		if(!is_eq_vector(prob_noisy.get_lb(), probs[i]->get_lb(), EPS) ||
		   !is_eq_vector(prob_noisy.get_ub(), probs[i]->get_ub(), EPS)){
			std::cout<<" bounds failed!"<<std::endl;
			return 1;
		}
		std::cout << " bounds passed. ";

		// Generate a test point (middle of the search space)
		decision_vector x(dim);
		for(unsigned int xi = 0; xi < x.size(); xi++){
			x[xi] = (probs[i]->get_lb()[xi] + probs[i]->get_ub()[xi]) / 2.0;
		}

		// Check that the sample mean of the noise is the same as specified by noise_mean
		fitness_vector f_noiseless = probs[i]->objfun(x);
		fitness_vector noise_sample_mean(f_noiseless.size(), 0.0);

		// Compute the sample mean
		// At the same time monitor if any of the noised points exceed the range
		// as defined by the noise params.
		for(int tr = 0; tr < num_trials; tr++){
			prob_noisy.set_seed(rng_seed_provider()); // Automatically resets the cache
			fitness_vector f_cur = prob_noisy.objfun(x);
			for(unsigned int fi = 0; fi < f_cur.size(); fi++){
				double diff_t = f_cur[fi] - f_noiseless[fi];
				noise_sample_mean[fi] +=  diff_t / (double)num_trials;	
				if(f_cur[fi] < f_noiseless[fi] + noise_lb){
					std::cout << " noise exceeds lower bound!" << std::endl;
					return 1;
				}
				if(f_cur[fi] > f_noiseless[fi] + noise_ub){
					std::cout << " noise exceeds upper bound!" << std::endl;
					return 1;
				}
			}
		}
		std::cout << "noise's upper & lower bounds passes. ";

		// Check if mean is sane under uniform noise
		double noise_mean = (noise_ub + noise_lb) / 2.0;
		for(unsigned int fi = 0; fi < noise_sample_mean.size(); fi++){
			if(fabs(noise_sample_mean[fi] - noise_mean) > tol){
				std::cout<< " mean of noise is strange!"<<std::endl;
				std::cout<< noise_sample_mean << ", expected = " << noise_mean << std::endl;
				return 1;
			}
		}
		std::cout << " mean of noise passed. " << std::endl;
	}

	// All good
	return 0;
}

int main()
{	
	int dimension = 10;
	std::vector<problem::base_ptr> probs;
	probs.push_back(problem::zdt(1, dimension).clone());
	probs.push_back(problem::ackley(dimension).clone());

	return test_noisy(probs, 0.0, 0.1, 5000, 0.01) ||
		   test_noisy(probs, 3.14, 0.1, 5000, 0.01) ||
		   test_noisy_uniform(probs, 0.0, 0.1, 5000, 0.01) ||
		   test_noisy_uniform(probs, -0.2, 0.2, 5000, 0.01);
}
