 /*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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
#include <iostream>
#include <iomanip>

#include "../src/pagmo.h"

using namespace pagmo;

/**
 * Generates a random square matrix with the number of vertices set to dimension.
 * @param dimension - the number of vertices in the matrix
 * @param verbose - prints the matrix to console if true
 * @return a square adjacency matrix
 */
std::vector<std::vector<double> > generate_random_matrix(int dimension, bool verbose = false) {
    std::default_random_engine rengine(time(NULL)); // seed software PRNG
    std::uniform_real_distribution<double> distr(0.1, 1); // range
    
    std::vector<std::vector<double> > random_2D(dimension, std::vector<double>(dimension, 0));

    if (verbose) std::cout << "Two dimensional matrix created:";
    
    for (int i = 0; i < dimension; ++i) {
        if (verbose) std::cout << std::endl;
        
        for (int j = 0; j < dimension; ++j) {
            if (i == j) 
                random_2D[i][j] = 0;
            else
                random_2D[i][j] = distr(rengine);

            if (verbose) std::cout << random_2D[i][j] << " \t "; 
        }
    }
    if (verbose) std::cout << std::endl ;
    
    return random_2D;
}

/*
 * This test creates a random vector<vector<double>> two dimensional vector,
 * then it instantiates a tsp object which loads the vector object.
 * The internal graph is then returned and used to instantiate a new problem.
 * The weights of the two tsp object matrices are then compared.
 * @param[in] repeat - the number of times to repeat the test
 * @param[in] l_bounds - the minimum random size of the square matrix
 * @param[in] u_bounds - the maximum random size of the square matrix
 * @param[in] verbose - print matrix and converted to console
 */
bool test_conversion(int repeat, int l_bounds, int u_bounds, bool verbose = false) {
    for (int i = 0; i < repeat; ++i) {
        // create random 2d vector and output it to console
        int no_vertices = rand() % u_bounds + l_bounds; // between l_b and u_b
        std::vector<std::vector<double> > original( generate_random_matrix(no_vertices, verbose) );

        // instantiate a tsp problem, vector constructor is called
        pagmo::problem::tsp prob(original);

        // output the graph structure, conversion done internally
        if (verbose)
            std::cout << prob.human_readable();

        // get the converted graph
        problem::tsp_graph graph = prob.get_graph();

        pagmo::problem::tsp new_prob(graph);

        // check equality
        if (prob.get_weights() != new_prob.get_weights() ) {
            std::cout << "vector2D to boost graph to vector2D conversion failed!\n";
            return true;
        }
    }
    return false;
}

int main()
{
    if (test_conversion(100, 10, 50, false)) return 1;
    
    // all iz well
    return 0;
}