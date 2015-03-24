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

// Test code for the shifted meta-problem

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include "../src/pagmo.h"
#include "test.h"

using namespace pagmo;

const double EPS = 10e-9;

int shifted_test1(){
	pagmo::problem::dtlz prob(2, 40);	
	decision_vector shift_vec(prob.get_dimension(), 6.0);
	pagmo::problem::shifted prob_shifted(prob, shift_vec);
	decision_vector p1_shifted_space(prob.get_dimension(), 6.2);
	decision_vector p1_original_space(prob.get_dimension(), 0.2);
	fitness_vector f_expected = prob.objfun(p1_original_space);
	fitness_vector f_shifted = prob_shifted.objfun(p1_shifted_space);
	constraint_vector c_expected = prob.compute_constraints(p1_original_space);
	constraint_vector c_shifted = prob_shifted.compute_constraints(p1_shifted_space);
	if(is_eq_vector(f_shifted, f_expected, EPS)){
		std::cout<<prob_shifted.get_name()<<" fitness passes, "<<"\t";
	}
	else{	
		std::cout<<prob_shifted.get_name()<<" fitness failed!"<<std::endl;
		PRINT_VEC(f_expected);
		PRINT_VEC(f_shifted);
		return 1;
	}
	if(is_eq_vector(c_shifted, c_expected, EPS)){
		std::cout<<" constraints passes."<<std::endl;
	}
	else{
		std::cout<<" constraints failed!"<<std::endl;
		PRINT_VEC(c_expected);
		PRINT_VEC(c_shifted);
		return 1;
	}	
	return 0;
}

int shifted_test2(){
	pagmo::problem::dtlz prob(1, 40);	
	decision_vector shift_vec(prob.get_dimension(), -6.0);
	pagmo::problem::shifted prob_shifted(prob, shift_vec);
	decision_vector p1_shifted_space(prob.get_dimension(), -5.8);
	decision_vector p1_original_space(prob.get_dimension(), 0.2);
	fitness_vector f_expected = prob.objfun(p1_original_space);
	fitness_vector f_shifted = prob_shifted.objfun(p1_shifted_space);
	constraint_vector c_expected = prob.compute_constraints(p1_original_space);
	constraint_vector c_shifted = prob_shifted.compute_constraints(p1_shifted_space);
	if(is_eq_vector(f_shifted, f_expected, EPS)){
		std::cout<<prob_shifted.get_name()<<" fitness passes, "<<"\t";
	}
	else{	
		std::cout<<prob_shifted.get_name()<<" fitness failed!"<<std::endl;
		PRINT_VEC(f_expected);
		PRINT_VEC(f_shifted);
		return 1;
	}
	if(is_eq_vector(c_shifted, c_expected, EPS)){
		std::cout<<" constraints passes."<<std::endl;
	}
	else{
		std::cout<<" constraints failed!"<<std::endl;
		PRINT_VEC(c_expected);
		PRINT_VEC(c_shifted);
		return 1;
	}	
	return 0;
}

int shifted_test3(){
	pagmo::problem::zdt prob(1, 40);	
	pagmo::problem::shifted prob_shifted(prob, -6); // Another interface
	decision_vector p1_shifted_space(prob.get_dimension(), -5.8);
	decision_vector p1_original_space(prob.get_dimension(), 0.2);
	fitness_vector f_expected = prob.objfun(p1_original_space);
	fitness_vector f_shifted = prob_shifted.objfun(p1_shifted_space);
	constraint_vector c_expected = prob.compute_constraints(p1_original_space);
	constraint_vector c_shifted = prob_shifted.compute_constraints(p1_shifted_space);
	if(is_eq_vector(f_shifted, f_expected, EPS)){
		std::cout<<prob_shifted.get_name()<<" fitness passes, "<<"\t";
	}
	else{	
		std::cout<<prob_shifted.get_name()<<" fitness failed!"<<std::endl;
		PRINT_VEC(f_expected);
		PRINT_VEC(f_shifted);
		return 1;
	}
	if(is_eq_vector(c_shifted, c_expected, EPS)){
		std::cout<<" constraints passes."<<std::endl;
	}
	else{
		std::cout<<" constraints failed!"<<std::endl;
		PRINT_VEC(c_expected);
		PRINT_VEC(c_shifted);
		return 1;
	}	
	return 0;
}

int shifted_test4(){
	pagmo::problem::ackley prob(40);	
	pagmo::problem::shifted prob_shifted(prob, -6); // Another interface
	decision_vector p1_shifted_space(prob.get_dimension(), -5.8);
	decision_vector p1_original_space(prob.get_dimension(), 0.2);
	fitness_vector f_expected = prob.objfun(p1_original_space);
	fitness_vector f_shifted = prob_shifted.objfun(p1_shifted_space);
	constraint_vector c_expected = prob.compute_constraints(p1_original_space);
	constraint_vector c_shifted = prob_shifted.compute_constraints(p1_shifted_space);
	if(is_eq_vector(f_shifted, f_expected, EPS)){
		std::cout<<prob_shifted.get_name()<<" fitness passes,"<<"\t";
	}
	else{	
		std::cout<<prob_shifted.get_name()<<" fitness failed!"<<std::endl;
		PRINT_VEC(f_expected);
		PRINT_VEC(f_shifted);
		return 1;
	}
	if(is_eq_vector(c_shifted, c_expected, EPS)){
		std::cout<<" constraints passes."<<std::endl;
	}
	else{
		std::cout<<" constraints failed!"<<std::endl;
		PRINT_VEC(c_expected);
		PRINT_VEC(c_shifted);
		return 1;
	}	
	return 0;
}

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
int test_shifted(
	const std::vector<problem::base_ptr> & probs,
	double d_from_center,
	double amount_shift)
{
	std::cout<<
	"Start batch testing with d_from_center = "<<
	d_from_center<<
	" and shift = "<<
	amount_shift<<
	std::endl;

	for(unsigned int i = 0; i < probs.size(); i++){
		pagmo::problem::shifted prob_shifted(*(probs[i]), amount_shift);
		decision_vector p_shifted_space = construct_test_point(prob_shifted.clone(), d_from_center);
		// Obtain the corresponding point in the original space
		decision_vector p_original_space(p_shifted_space.size(), 0);
		for(unsigned int k = 0; k < p_original_space.size(); k++){
			p_original_space[k] = p_shifted_space[k] - amount_shift;
		}
		decision_vector f_shifted = prob_shifted.objfun(p_shifted_space);
		decision_vector f_original = probs[i]->objfun(p_original_space);
		decision_vector c_shifted = prob_shifted.compute_constraints(p_shifted_space);
		decision_vector c_original = probs[i]->compute_constraints(p_original_space); 
		if(is_eq_vector(f_shifted, f_original, EPS)){
			std::cout << std::setw(40) << prob_shifted.get_name() << " fitness passes, ";
		}
		else{	
			std::cout << prob_shifted.get_name() << " fitness failed!"<<std::endl;
			PRINT_VEC(f_original);
			PRINT_VEC(f_shifted);
			return 1;
		}
		if(is_eq_vector(c_shifted, c_original, EPS)){
			std::cout << " constraints passes." << std::endl;
		}
		else{
			std::cout <<" constraints failed!" <<std::endl;
			PRINT_VEC(c_original);
			PRINT_VEC(c_shifted);
			return 1;
		}	
	}
	return 0;
}

int main()
{		
	int dimension = 40;
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
	return shifted_test1() ||
		   shifted_test2() ||
		   shifted_test3() ||
		   shifted_test4() ||
		   test_shifted(probs, 0.2, 999) ||
		   test_shifted(probs, 0.2, -999) ||
		   test_shifted(probs, -0.2, 0) ||
		   test_shifted(probs, 0, 0) ||
		   test_shifted(probs, -0.99, 3.14) ||
		   test_shifted(probs, 0.99, -3.14);
}
