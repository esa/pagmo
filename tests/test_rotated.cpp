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

// Test code for the rotated meta-problem

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include "../src/pagmo.h"
#include "../src/Eigen/Dense"
#include "test.h"

using namespace pagmo;

const double EPS = 10e-9;

// Construct a test point in the new transformed space
decision_vector construct_test_point(const problem::base_ptr &prob, double d_from_center)
{
	assert(d_from_center <= 1.0);
	decision_vector lb = prob->get_lb();
	decision_vector ub = prob->get_ub();
	decision_vector test_point(prob->get_dimension(), 0);
	for(unsigned int i = 0; i < prob->get_dimension(); i++){
		if(abs(ub[i] - lb[i]) < EPS) continue;
		double middle_t = (ub[i] + lb[i]) / 2.0;
		test_point[i] = middle_t + d_from_center * (ub[i] - middle_t);
	}
	return test_point;
}

// Run the batch test by constructing meta-problems against probs
// d_from_center: Percentage of "deviation" away from the middle of the bounds
// d_from_center = -1 or 1 means at the boundary points
// scenario: [0 | 1]; 0: Identity matrix; 1: Random orthogonal matrix
int test_rotated(
	const std::vector<problem::base_ptr> & probs,
	double d_from_center,
	int scenario)
{
	std::cout <<
	"Start batch testing with d_from_center = " <<
	d_from_center << " with " << 
	(scenario==0?"identity matrix.":"random orthogonal matrix.") <<
	std::endl;

	for(unsigned int i = 0; i < probs.size(); i++){

		int dim = probs[i]->get_dimension();

		Eigen::MatrixXd Rot;
		if(scenario==0){
			Rot = Eigen::MatrixXd::Identity(dim, dim);
		}else{
			Rot = Eigen::MatrixXd::Random(dim, dim).householderQr().householderQ();
		}

		pagmo::problem::rotated prob_rotated(*(probs[i]), Rot);

		std::cout<< std::setw(40) << prob_rotated.get_name();

		// Check that bounds are normalized
		for(unsigned k = 0; k < prob_rotated.get_lb().size(); k++){
			if(fabs(prob_rotated.get_lb()[k] - (-sqrt(2.0))) > EPS ||
			   fabs(prob_rotated.get_ub()[k] - sqrt(2.0)) > EPS){
				std::cout<<" bounds failed!"<<std::endl;
				return 1;
			}
		}
		std::cout << " bounds passed. ";

		decision_vector p_rotated_space = construct_test_point(prob_rotated.clone(), d_from_center);

		// Obtain the corresponding point in the original space
		decision_vector p_original_space(p_rotated_space.size(), 0);
			
		Eigen::VectorXd p_vec_rotated_space(p_rotated_space.size());
		Eigen::VectorXd p_vec_original_space;
		for(unsigned int k = 0; k < p_rotated_space.size(); k++){
			p_vec_rotated_space(k) = p_rotated_space[k];
		}
		// De-rotation
		p_vec_original_space = Rot.transpose() * p_vec_rotated_space;
		// De-normalize
		for(unsigned int k = 0; k < p_original_space.size(); k++){
			double ub_t = probs[i]->get_ub()[k];
			double lb_t = probs[i]->get_lb()[k];
			p_original_space[k] = p_vec_original_space(k)*((ub_t-lb_t)/2.0) + (ub_t+lb_t)/2.0;
			// Projection
			p_original_space[k] = std::min(ub_t, p_original_space[k]);	
			p_original_space[k] = std::max(lb_t, p_original_space[k]);
		}

		fitness_vector f_rotated = prob_rotated.objfun(p_rotated_space);
		fitness_vector f_original = probs[i]->objfun(p_original_space);
		constraint_vector c_rotated = prob_rotated.compute_constraints(p_rotated_space);
		constraint_vector c_original = probs[i]->compute_constraints(p_original_space); 

		if(is_eq_vector(f_rotated, f_original, EPS)){
			std::cout << " fitness passes, ";
		}
		else{	
			std::cout << " fitness failed!"<<std::endl;
			std::cout << "original fitness: " << f_original << std::endl;
			std::cout << "new fitness: " << f_rotated << std::endl;
			return 1;
		}
		if(is_eq_vector(c_rotated, c_original, EPS)){
			std::cout << " constraints passes." << std::endl;
		}
		else{
			std::cout <<" constraints failed!" <<std::endl;
			std::cout << "original constraints: " << c_original << std::endl;
			std::cout << "new constraints: " << c_rotated << std::endl;
			return 1;
		}	
	}
	return 0;
}

int main()
{		
	int dimension = 10;
	std::vector<problem::base_ptr> probs;

	probs.push_back(problem::zdt(1,dimension).clone());
	probs.push_back(problem::zdt(2,dimension).clone());
	probs.push_back(problem::zdt(3,dimension).clone());
	probs.push_back(problem::zdt(4,dimension).clone());
	probs.push_back(problem::zdt(6,dimension).clone());

	for (int i = 1;i <= 7; i++) {
	    probs.push_back(problem::dtlz(i,dimension).clone());
	}
	
	probs.push_back(problem::ackley(dimension).clone());
	probs.push_back(problem::rastrigin(dimension).clone());

	return test_rotated(probs, 0.0, 0) ||
		   test_rotated(probs, -0.1, 0) ||
		   test_rotated(probs, -0.2, 0) ||
		   test_rotated(probs, 0.1, 0) ||
		   test_rotated(probs, 0.2, 0) ||
		   test_rotated(probs, 0.0, 1) ||
		   test_rotated(probs, -0.1, 1) ||
		   test_rotated(probs, -0.2, 1) ||
		   test_rotated(probs, 0.1, 1) ||
		   test_rotated(probs, 0.2, 1) ||
		   test_rotated(probs, -0.5, 1) ||
		   test_rotated(probs, 0.5, 1);
}
