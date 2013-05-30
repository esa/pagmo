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
#include "../src/keplerian_toolbox/planet_ss.h"
#include "../src/keplerian_toolbox/epoch.h"

#include "../src/Eigen/Dense"

//--------------------------------------------------------------------------------
// static data needed to test the non-default constructor in some of the problems.
#ifdef PAGMO_ENABLE_KEP_TOOLBOX 
//mga_1dsm
const std::vector<kep_toolbox::planet_ptr> construct_sequence() {
	std::vector<kep_toolbox::planet_ptr> retval;
	retval.push_back(kep_toolbox::planet_ss("earth").clone());
	retval.push_back(kep_toolbox::planet_ss("earth").clone());
	retval.push_back(kep_toolbox::planet_ss("earth").clone());
	return retval;
};
#endif

//knapsack
static std::vector<double> a(11,30), b(11,10);
static double c = 15;

//laplace
static const int default_sequence[5] = {3,2,2,1,5};
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
	unsigned int dimension = 24;

	// create a containers of pagmo::problems
	std::vector<problem::base_ptr> probs;

	// fill it up with problems
	probs.push_back(problem::ackley(dimension).clone());
	probs.push_back(problem::rosenbrock(dimension).clone());
	probs.push_back(problem::branin().clone());
	probs.push_back(problem::dejong(dimension).clone());
	probs.push_back(problem::fon().clone());
	probs.push_back(problem::golomb_ruler(10,20).clone());
	probs.push_back(problem::griewank(dimension).clone());
	probs.push_back(problem::himmelblau().clone());
	probs.push_back(problem::string_match("e dai dai dai.....portiamolo a casa!!").clone());
	probs.push_back(problem::inventory(7,8,1234).clone());
	probs.push_back(problem::knapsack(a,b,c).clone());
	probs.push_back(problem::kur(dimension).clone());
	probs.push_back(problem::lennard_jones(dimension).clone());
	probs.push_back(problem::levy5(dimension).clone());
	probs.push_back(problem::luksan_vlcek_1(dimension).clone());
	probs.push_back(problem::luksan_vlcek_2(dimension).clone());
	probs.push_back(problem::luksan_vlcek_3(dimension).clone());
	probs.push_back(problem::michalewicz(dimension).clone());
	probs.push_back(problem::pol().clone());
	probs.push_back(problem::rastrigin(dimension).clone());
	probs.push_back(problem::sch().clone());
	probs.push_back(problem::schwefel(dimension).clone());
	probs.push_back(problem::snopt_toyprob().clone());
	probs.push_back(problem::zdt1(dimension).clone());
	probs.push_back(problem::zdt2(dimension).clone());
	probs.push_back(problem::zdt3(dimension).clone());
	probs.push_back(problem::zdt4(dimension).clone());
	probs.push_back(problem::zdt5(dimension).clone());
	probs.push_back(problem::zdt6(dimension).clone());
	probs.push_back(problem::dtlz1(dimension).clone());
	probs.push_back(problem::dtlz2(dimension).clone());
	probs.push_back(problem::dtlz3(dimension).clone());
	probs.push_back(problem::dtlz4(dimension).clone());
	probs.push_back(problem::dtlz5(dimension).clone());
	probs.push_back(problem::dtlz6(dimension).clone());
	probs.push_back(problem::dtlz7(dimension).clone());
	probs.push_back(problem::tsp().clone()); //TODO: define the tsp using a non-default weight-matrix

    //----- Test CEC2006 -----//
    for(int i=1; i<=24; i++){
        probs.push_back(problem::cec2006(i).clone());
    }

    //----- Test CEC2009 - UF set-----//
    for(int i=1; i<=10; i++){
		probs.push_back(problem::cec2009(i, dimension, false).clone());
    }
    //----- Test CEC2009 - CF set-----//
    for(int i=1; i<=10; i++){
		probs.push_back(problem::cec2009(i, dimension, true).clone());
	}

	//----- Test meta-problems -----//
	problem::zdt1 zdt1_before_transform1(dimension);
	//----- shifted -----//
	probs.push_back(problem::shifted(zdt1_before_transform1).clone());
	//----- rotated -----//
	probs.push_back(problem::rotated(zdt1_before_transform1).clone());

#ifdef PAGMO_ENABLE_KEP_TOOLBOX
	probs.push_back(problem::cassini_1(2).clone());
	probs.push_back(problem::cassini_2().clone());
	probs.push_back(problem::gtoc_1().clone());
	probs.push_back(problem::messenger().clone());
	probs.push_back(problem::rosetta().clone());
	probs.push_back(problem::messenger_full().clone());
	probs.push_back(problem::tandem(3,10).clone());
	probs.push_back(problem::laplace(std::vector<int>(default_sequence,default_sequence + 5)).clone());
	probs.push_back(problem::mga_1dsm_alpha(construct_sequence()).clone());
	probs.push_back(problem::mga_1dsm_tof(construct_sequence()).clone());
#endif	

    // initialize the best solution and retrieve the decision vector
    for (size_t i=0; i<probs.size(); ++i) {
        std::cout << std::endl << std::setw(40) << probs[i]->get_name() << std::endl;

        if(probs[i]->get_best_known_x_vector().empty())
            std::cout << "Best decision vector is not implemented." << std::endl;
        else {
            const std::vector<decision_vector> &x_best_known_vector = probs[i]->get_best_known_x_vector();

            for(int j=0; j<x_best_known_vector.size(); j++) {
                const decision_vector &x_best_known = x_best_known_vector.at(j);

                fitness_vector f_computed = probs[i]->objfun(x_best_known);
                const fitness_vector &f_best_known = probs[i]->get_best_known_f_vector().at(j);

                if(is_eq(f_computed, f_best_known, EPS)){
                    std::cout << " fitness passes, ";
                }
                else{
                    std::cout << " fitness failed!"<<std::endl;
                    return 1;
                }

                constraint_vector c_computed = probs[i]->compute_constraints(x_best_known);
                const constraint_vector &c_best_known = probs[i]->get_best_known_c_vector().at(j);

                if(c_best_known.empty()) {
                    std::cout << "Best constraint vector is not implemented." << std::endl;
                }
                else {
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

	std::cout << std::endl;
	return 0;
}
