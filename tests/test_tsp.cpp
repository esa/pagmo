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

/*
 * This test file creates a random vector<vector<double>> two dimensional vector,
 * then it instantiates a tsp object which loads the vector object.
 * It then internally converts the vector to a boost graph,
 * which is then returned and checked against the initial 2D vector.
 */

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

int main()
{
    // create random 2d vector and output it to console
    int no_vertices = rand() % 10 + 1; // between 1 and 10
    problem::vector2D<double> mat( generate_random_vector2D(no_vertices, true) );
    
    // instantiate a tsp problem, vector constructor is called
    pagmo::problem::tsp tsprob(mat);
    
    // output the graph structure, conversion done internally
    std::cout << tsprob.human_readable();
    
    // get the converted graph
    problem::tsp_graph g = tsprob.get_graph();
    
    // convert back to vector2D and check equality
    if (mat != pagmo::problem::base_tsp::graph_to_vector2D( g )) {
        std::cout << "vector2D to boost graph to vector2D conversion failed!\n";
        return 1;
    }
    
    return 0;
}