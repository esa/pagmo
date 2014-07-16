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

problem::vector2D<double> generate_random_vector2D(int dimension, bool verbose = false) {
    std::default_random_engine rengine(time(NULL)); // seed software PRNG
    std::uniform_real_distribution<> distr(0, 1); // range
    
    problem::vector2D<double> random_2D(dimension, std::vector<double>(dimension, 0));

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
 * It then internally converts the vector to a boost graph,
 * which is then returned and checked against the initial 2D vector.
 * @param[in] repeat - the number of times to repeat the test
 * @param[in] l_bounds - the minimum random size of the square matrix
 * @param[in] u_bounds - the maximum random size of the square matrix
 * @param[in] verbose - print matrix and converted to console
 */
bool test_conversion(int repeat, int l_bounds, int u_bounds, bool verbose = false) {
    for (int i = 0; i < repeat; ++i) {
        // create random 2d vector and output it to console
        int no_vertices = rand() % u_bounds + l_bounds; // between l_b and u_b
        problem::vector2D<double> original( generate_random_vector2D(no_vertices, verbose) );

        // instantiate a tsp problem, vector constructor is called
        pagmo::problem::tsp tsprob(original);

        // output the graph structure, conversion done internally
        if (verbose)
            std::cout << tsprob.human_readable();

        // get the converted graph
        problem::tsp_graph graph = tsprob.get_graph();

        // convert back to vector2D
        problem::vector2D<double> converted(no_vertices, std::vector<double>(no_vertices, 0));
        pagmo::problem::base_tsp::convert_graph_to_vector2D(graph, converted);

        // check equality
        if (original != converted) {
            std::cout << "vector2D to boost graph to vector2D conversion failed!\n";
            return true;
        }
    }
    return false;
}

/**
 * This tests the compute_idx function, which returns the index in the
 * one dimensional vector of concatenated rows of a matrix (vector2D)
 * For testing we instantiate a matrix with 4 rows and then compute the
 * sum for the columns and rows, skipping elements from the main diagonal.
 * 
 * Create test matrix with values from 0 to n^2-1
 * matrix = 
 * ______________
 * | 0  1  2  3 |->  1 +  2 +  3 = 6
 * | 4  5  6  7 |->  4 +  6 +  7 = 17
 * | 8  9 10 11 |->  8 +  9 + 11 = 28
 * |12 13 14 15 |-> 12 + 13 + 14 = 39
 * |_|__|__|__|_|
 *   |--|--|--|--->  4 +  8 + 12 = 24
 *      |--|--|--->  1 +  9 + 13 = 23
 *         |--|--->  2 +  6 + 14 = 22
 *            |--->  3 +  7 + 11 = 21
 * 
 * then v = {1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14}, where |v| = n*(n-1) = 12
 * 
 * and c = {6, 17, 28, 39, 24, 23, 22, 21, 0, 0, 0, 0, 0, 0}, 
 *      where |c!=0| = 2*n = 8
 * 
 * Inequalities: length = (n-1)(n-2)
 * 
 * c[in] = ui - uj + n*x_{i,j} - n
 * 
 *   |uj 1  2  3  4
 * __|_____________
 * ui|
 * 1 |   -  -  -  -
 * 2 |   -  -  6  7
 * 3 |   -  9  - 11
 * 4 |   - 13 14  -
 * 
 * 2 - 3 + 4 *  6 - 4 = 19
 * 2 - 4 + 4 *  7 - 4 = 22
 * 3 - 2 + 4 *  9 - 4 = 33
 * 3 - 4 + 4 * 11 - 4 = 39
 * 4 - 2 + 4 * 13 - 4 = 50
 * 4 - 3 + 4 * 14 - 4 = 53
 * 
 * now c = {6, 17, 28, 39, 24, 23, 22, 21, 18, 22, 33, 39, 50, 53},
 *      where |c| = 2*n + (n-1)*(n-2) = 14
 * 
 * @param[in] verbose - prints matrix, indexes and resulting c
 */
bool test_compute_idx(bool verbose = false) {
    int n = 4; // the number of vertices for the square matrix
    int ceq = 2*n; // the number of equality constraints
    problem::vector2D<int> matrix(n, std::vector<int>(n, 0));
    std::vector<int> x(n*(n-1), 0);
    std::vector<int> c(2*n + (n-1)*(n-2), 0);
    // create check, must be finally equal to c
    std::vector<int> check = {6, 17, 28, 39, 24, 23, 22, 21, 19, 22, 33, 39, 50, 53};
    
    int k = 0;
    // create matrix and v
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = i*n + j;
            if (i!=j) 
                x[k++] = matrix[i][j];
        }
    }
    
    // compute row and col sums
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if(i==j) continue; // ignoring main diagonal

            int idx_row = problem::tsp::compute_idx(i, j, n);
            int idx_col = problem::tsp::compute_idx(j, i, n);

            // equalities
            c[i  ] += x[idx_row];
            c[i+n] += x[idx_col];
            
            // inequalities ( ignoring first row & column )
            if(i != 0 && j != 0)
                c[ceq++] = (i+1) - (j+1) + n * x[idx_row] - n;
            
            if (verbose)
                std::cout << i << " - " << j << " = " << matrix[i][j] 
                            << " \t idx_row: " << idx_row
                            << " \t v[idx_row]: " << x[idx_row]    
                            << " \t idx_col: " << idx_col
                            << " \t v[idx_col]: " << x[idx_col]
                            << std::endl;
        }
    }
    
    if (verbose)
        for (int i = 0; i < (int)c.size(); ++i)
            std::cout << c[i] << " ";
    
    // check equality
    if (check != c) {
        std::cout << "compute_idx function applied to row & column sums failed!\n";
        return true;
    }
    return false;
}

int main()
{
    if (!test_conversion(100, 10, 20, true)) return 1;
    if (!test_compute_idx(true)) return 1;
    
    return 0;
}