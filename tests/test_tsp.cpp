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
#include <iostream>
#include <iomanip>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"

#include "../src/problem/tsp.h"
#include "../src/population.h"

using namespace pagmo;

/**
 * Generates a random square matrix dimension x dimension.
 * @param dimension - the matrix dimension
 * @return a square adjacency matrix
 */
std::vector<std::vector<double> > generate_random_matrix(int dimension, boost::lagged_fibonacci607 rng) {

    boost::uniform_real<double> uniform(0.0,1.0);
    boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > distr(rng,uniform);
   
    std::vector<std::vector<double> > retval(dimension, std::vector<double>(dimension, 0));
  
    for (int i = 0; i < dimension; ++i) {        
        for (int j = 0; j < dimension; ++j) {
            if (i == j) 
                retval[i][j] = 0;
            else
                retval[i][j] = distr();
        }
    }
    return retval;
}

/*
 * This test creates three tsp problems (one per encoding type) and
 * checks that objective function and feasibility is invariant across problem
 * representation
 *
 * @param[in] repeat - the number of times to repeat the test
 * @param[in] l_bounds - the minimum random size of the square matrix
 * @param[in] u_bounds - the maximum random size of the square matrix
 * @param[in] verbose - print matrix and converted to console
 */
bool test_encodings_equivalence(int repeat, int l_bounds, int u_bounds, boost::lagged_fibonacci607 rng) 
{
    for (int i = 0; i < repeat; ++i) {
        // create random 2d vector and output it to console
        boost::uniform_int<int> uniform(l_bounds,u_bounds);
        boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_int<int> > distr(rng,uniform);
        int n_cities = distr(); // between l_b and u_b
        std::vector<std::vector<double> > weights( generate_random_matrix(n_cities,rng) );

        // instantiate a tsp problem for each of the available encodings
        pagmo::problem::tsp prob_full(weights, pagmo::problem::tsp::FULL);
        pagmo::problem::tsp prob_rk(weights, pagmo::problem::tsp::RANDOMKEYS);
        pagmo::problem::tsp prob_cities(weights, pagmo::problem::tsp::CITIES);

        pagmo::decision_vector tour_rk = population(prob_rk,1).get_individual(0).cur_x;
        pagmo::decision_vector tour_cities = prob_rk.randomkeys2cities(tour_rk);
        pagmo::decision_vector tour_full = prob_full.cities2full(tour_cities);

        pagmo::fitness_vector f_rk = prob_rk.objfun(tour_rk);
        pagmo::fitness_vector f_cities = prob_cities.objfun(tour_cities);
        pagmo::fitness_vector f_full = prob_full.objfun(tour_full);

        // check equality
        if ( (f_rk!=f_cities) || (f_rk!=f_full) ) {
            std::cout << "fitness is different across encodings\n";
            return true;
        }
        if ( (!prob_full.feasibility_x(tour_full)) || (!prob_cities.feasibility_x(tour_cities)) || (!prob_rk.feasibility_x(tour_rk)) ) 
        {
            std::cout << "feasibility is different across encodings\n";
            return true;
        }
    }
    return false;
}


/*
 * This test creates three tsp problems (one per encoding type) and
 * checks that objective function and feasibility is invariant across problem
 * representation
 *
 * @param[in] repeat - the number of times to repeat the test
 * @param[in] l_bounds - the minimum random size of the square matrix
 * @param[in] u_bounds - the maximum random size of the square matrix
 * @param[in] verbose - print matrix and converted to console
 */
bool test_encoding_transformations(int repeat, boost::lagged_fibonacci607 rng) 
{
    for (int i = 0; i < repeat; ++i) {
        // create random 2d vector and output it to console
        boost::uniform_int<int> uniform(3,50);
        boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_int<int> > distr(rng,uniform);
        int n_cities = distr(); // between l_b and u_b
        std::vector<std::vector<double> > weights( generate_random_matrix(n_cities,rng) );

        pagmo::problem::tsp prob_rk(weights, pagmo::problem::tsp::RANDOMKEYS);
        pagmo::decision_vector tour_rk = population(prob_rk,1).get_individual(0).cur_x;

        pagmo::fitness_vector f1 = prob_rk.objfun( tour_rk );
        pagmo::fitness_vector f2 = prob_rk.objfun( prob_rk.cities2randomkeys(prob_rk.full2cities(prob_rk.cities2full(prob_rk.randomkeys2cities(tour_rk) ) ), tour_rk ) );
        if (f1!=f2) 
        {
            return true;
        }
    }
    return false;
}

int main()
{
    boost::lagged_fibonacci607 rng;
    std::cout << "Testing Encoding Equivalence: ";
    if (test_encodings_equivalence(100, 3, 50, rng)) return 1;
    std::cout << "SUCCESS" << std::endl;
    std::cout << "Testing Encoding Transformations: ";
    if (test_encoding_transformations(100,rng)) return 1;
    std::cout << "SUCCESS" << std::endl;
    
    // all iz well
    return 0;
}