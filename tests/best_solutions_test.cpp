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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "../src/pagmo.h"

//--------------------------------------------------------------------------------

const double EPS = 1e-9;

///The idea of this unit test is to go through all pagmo::problems, initialize the best solution,
///retrieve the decision vector 

using namespace pagmo;

bool is_eq(const fitness_vector & f1, const fitness_vector & f2, double eps){
    if(f1.size() != f2.size()) return false;
    for(unsigned int i = 0; i < f1.size(); i++){
        if(fabs(f1[i]-f2[i])>eps) return false;
    }
    return true;
}

int main()
{
    //----- Test CEC2006 -----//
    std::vector<problem::cec2006*> probs;

    //----- Test CEC2006 -----//
    for(int i=1; i<=24; i++){
        probs.push_back(new problem::cec2006(i));
    }

    // initialize the best solution and retrieve the decision vector
    for (size_t i=0; i<probs.size(); ++i) {
        std::cout << std::endl << std::setw(40) << probs[i]->get_name() << std::endl;

        if(boost::math::isnan(((probs[i]->get_best_known_x_vector())[0])))
            std::cout << "Best decision vector is not implemented." << std::endl;
        else {
            const std::vector<decision_vector> &x_best_known_vector = probs[i]->get_best_known_x_vector();

            for(int j=0; j<x_best_known_vector.size(); j++) {
                const decision_vector &x_best_known = x_best_known_vector.at(j);

                fitness_vector f_computed = probs[i]->objfun(x_best_known);
                constraint_vector c_computed = probs[i]->compute_constraints(x_best_known);

                fitness_vector f_best_known = probs[i]->get_best_known_f_vector().at(j);
                constraint_vector c_best_known = probs[i]->get_best_known_c_vector().at(j);

                if(is_eq(f_computed, f_best_known, EPS)){
                    std::cout << " fitness passes, ";
                }
                else{
                    std::cout << " fitness failed!"<<std::endl;
                    return 1;
                }

                bool check_constraints = true;
                for(int j=0; j<c_computed.size(); j++)
                {
                    if(boost::math::isnan(c_computed.at(j))) {
                        std::cout << "Best constraint vector is not implemented." << std::endl;
                        check_constraints = false;
                        break;
                    }
                }

                if(check_constraints) {
                    if(is_eq(c_computed, c_best_known, EPS)){
                        std::cout << " constraints passes.";
                    }
                    else{
                        std::cout << " constraints failed!"<<std::endl;
                        return 1;
                    }
                }
            }
        }
    }

    //----- Delete test CEC2006 -----//
    for(int i=1; i<=24; i++){
        delete probs[i];
    }

	std::cout << std::endl;
	return 0;
}
