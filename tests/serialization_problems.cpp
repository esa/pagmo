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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#ifdef PAGMO_ENABLE_KEP_TOOLBOX
	#include <keplerian_toolbox/planet/jpl_low_precision.h>
	#include <keplerian_toolbox/epoch.h>
#endif

#include "../src/pagmo.h"

#include "../src/Eigen/Dense"
#include <boost/shared_ptr.hpp>
//-------------------------------------------------------------------------------
// static data needed to test the non-default constructor in some of the problems.
#ifdef PAGMO_ENABLE_KEP_TOOLBOX
//mga_1dsm
const std::vector<kep_toolbox::planet::planet_ptr> construct_sequence() {
	std::vector<kep_toolbox::planet::planet_ptr> retval;
	retval.push_back(kep_toolbox::planet::jpl_lp("earth").clone());
	retval.push_back(kep_toolbox::planet::jpl_lp("earth").clone());
	retval.push_back(kep_toolbox::planet::jpl_lp("earth").clone());
	return retval;
}
#endif


//laplace
static const int default_sequence[5] = {3,2,2,1,5};
//--------------------------------------------------------------------------------

///The idea of this unit test is to serialize all pagmo::problems, deserialize them and check that
///the objective function and the constraint implementation return the same in the original and in the deserialized object

using namespace pagmo;
int main()
{
	unsigned int dimension = 24;

	// create two containers of pagmo::problems
	std::vector<problem::base_ptr> probs;
	std::vector<problem::base_ptr> probs_new;

	// fill it up with problems
	probs.push_back(problem::ackley(dimension).clone());
	probs_new.push_back(problem::ackley().clone());
	probs.push_back(problem::rosenbrock(dimension).clone());
	probs_new.push_back(problem::rosenbrock().clone());
	probs.push_back(problem::branin().clone());
	probs_new.push_back(problem::branin().clone());
	probs.push_back(problem::dejong(dimension).clone());
	probs_new.push_back(problem::dejong().clone());
	probs.push_back(problem::fon().clone());
	probs_new.push_back(problem::fon().clone());
	probs.push_back(problem::golomb_ruler(10,20).clone());
	probs_new.push_back(problem::golomb_ruler().clone());
	probs.push_back(problem::griewank(dimension).clone());
	probs_new.push_back(problem::griewank().clone());
	probs.push_back(problem::himmelblau().clone());
	probs_new.push_back(problem::himmelblau().clone());
	probs.push_back(problem::string_match("e dai dai dai.....portiamolo a casa!!").clone());
	probs_new.push_back(problem::string_match().clone());
	probs.push_back(problem::inventory(7,8,1234).clone());
	probs_new.push_back(problem::inventory().clone());
	probs.push_back(problem::kur(dimension).clone());
	probs_new.push_back(problem::kur().clone());
	probs.push_back(problem::lennard_jones(dimension).clone());
	probs_new.push_back(problem::lennard_jones().clone());
	probs.push_back(problem::lavor_maculan(dimension).clone());
	probs_new.push_back(problem::lavor_maculan().clone());
	probs.push_back(problem::levy5(dimension).clone());
	probs_new.push_back(problem::levy5().clone());
	probs.push_back(problem::luksan_vlcek_1(dimension).clone());
	probs_new.push_back(problem::luksan_vlcek_1().clone());
	probs.push_back(problem::luksan_vlcek_2(dimension).clone());
	probs_new.push_back(problem::luksan_vlcek_2().clone());
	probs.push_back(problem::luksan_vlcek_3(dimension).clone());
	probs_new.push_back(problem::luksan_vlcek_3().clone());
	probs.push_back(problem::michalewicz(dimension).clone());
	probs_new.push_back(problem::michalewicz().clone());
	probs.push_back(problem::pol().clone());
	probs_new.push_back(problem::pol().clone());
	probs.push_back(problem::rastrigin(dimension).clone());
	probs_new.push_back(problem::rastrigin().clone());
	probs.push_back(problem::sch().clone());
	probs_new.push_back(problem::sch().clone());
	probs.push_back(problem::schwefel(dimension).clone());
	probs_new.push_back(problem::schwefel().clone());
	probs.push_back(problem::snopt_toyprob().clone());
	probs_new.push_back(problem::snopt_toyprob().clone());
	
    //----- Test TSP -----//
    probs_new.push_back(problem::tsp().clone());

	//----- Test ZDT -----//
	for(int i = 1; i <= 6; i++) {
	probs.push_back(problem::zdt(i, dimension).clone());
	probs_new.push_back(problem::zdt(i%6 + 1).clone());
	}

	//----- Test DTLZ -----//
	for(int i = 1; i <= 7; i++) {
		probs.push_back(problem::dtlz(i, dimension, 10).clone());
		probs_new.push_back(problem::dtlz(i%7 + 1).clone());
	}

	//----- Test CEC2006 -----//
	for(int i=1; i<=24; i++){
		probs.push_back(problem::cec2006(i).clone());
		probs_new.push_back(problem::cec2006(i%24 + 1).clone());
	}

	//----- Test CEC2009 - UF set-----//
	for(int i=1; i<=10; i++){
		probs.push_back(problem::cec2009(i, dimension, false).clone());
		probs_new.push_back(problem::cec2009(i%10 + 1, 11, false).clone());
	}
	//----- Test CEC2009 - CF set-----//
	for(int i=1; i<=10; i++){
		probs.push_back(problem::cec2009(i, dimension, true).clone());
		probs_new.push_back(problem::cec2009(1%10 + 1, 17, true).clone());
	}

	//----- Test meta-problems -----//
	problem::zdt zdt1_before_transform1(1, dimension);
	//----- shifted -----//
	probs.push_back(problem::shifted(zdt1_before_transform1).clone());
	probs_new.push_back(problem::shifted(zdt1_before_transform1).clone());
	//----- rotated -----//
	probs.push_back(problem::rotated(zdt1_before_transform1).clone());
	probs_new.push_back(problem::rotated(zdt1_before_transform1).clone()); //Will have a different random rotation matrix
	//----- noisy -----//
	probs.push_back(problem::noisy(zdt1_before_transform1,0,0,1.0,
				problem::noisy::NORMAL).clone());
	probs_new.push_back(problem::noisy(zdt1_before_transform1,111,1.0,4.5,
					problem::noisy::UNIFORM).clone());
	//----- robust ----- //
	probs.push_back(problem::robust(zdt1_before_transform1, 10, 0.1, 123).clone());
	probs_new.push_back(problem::robust(zdt1_before_transform1, 1, 1.23, 456).clone());

	//----- Test constraints handling meta-problems -----//
	problem::cec2006 cec2006_before_cstrs_handling(7);
	probs.push_back(problem::cstrs_self_adaptive(cec2006_before_cstrs_handling).clone());
	probs_new.push_back(problem::cstrs_self_adaptive(cec2006_before_cstrs_handling).clone());

	probs.push_back(problem::death_penalty(cec2006_before_cstrs_handling,problem::death_penalty::KURI).clone());
	probs_new.push_back(problem::death_penalty(cec2006_before_cstrs_handling,problem::death_penalty::SIMPLE).clone());

#ifdef PAGMO_ENABLE_KEP_TOOLBOX
	probs.push_back(problem::cassini_1(2).clone());
	probs_new.push_back(problem::cassini_1().clone());
	probs.push_back(problem::cassini_2().clone());
	probs_new.push_back(problem::cassini_2().clone());
	probs.push_back(problem::gtoc_1().clone());
	probs_new.push_back(problem::gtoc_1().clone());
	probs.push_back(problem::messenger().clone());
	probs_new.push_back(problem::messenger().clone());
	probs.push_back(problem::rosetta().clone());
	probs_new.push_back(problem::rosetta().clone());
	probs.push_back(problem::messenger_full().clone());
	probs_new.push_back(problem::messenger_full().clone());
	probs.push_back(problem::tandem(3,10).clone());
	probs_new.push_back(problem::tandem().clone());
	probs.push_back(problem::laplace(std::vector<int>(default_sequence,default_sequence + 5)).clone());
	probs_new.push_back(problem::laplace().clone());
	probs.push_back(problem::mga_1dsm_alpha(construct_sequence()).clone());
	probs_new.push_back(problem::mga_1dsm_alpha().clone());
	probs.push_back(problem::mga_1dsm_tof(construct_sequence()).clone());
	probs_new.push_back(problem::mga_1dsm_tof().clone());
#endif

	//serialize probs and deserialize into probs_new checking they are then identical
	for (size_t i=0; i< probs.size(); ++i) {
		{
		// create and open a character archive for output
		std::ofstream ofs("test.ar");
		// save data to archive
		boost::archive::text_oarchive oa(ofs);
		// write class instance to archive
		oa & probs[i];
		// archive and stream closed when destructors are called
		}

		{
		// create and open an archive for input
		std::ifstream ifs("test.ar");
		boost::archive::text_iarchive ia(ifs);
		// read class state from archive
		ia & probs_new[i];
		// archive and stream closed when destructors are called
		}

		{
		std::cout << std::endl << std::setw(40) << probs[i]->get_name()<<std::flush;
		decision_vector x(probs[i]->get_dimension(),0);
		fitness_vector f1(probs[i]->get_f_dimension(),0), f2(probs[i]->get_f_dimension(),1);
		constraint_vector c1(probs[i]->get_c_dimension(),0), c2(probs[i]->get_c_dimension(),1);
		population pop(*probs[i],1);
		x = pop.champion().x;
		probs[i]->objfun(f1,x);
		probs_new[i]->objfun(f2,x);
		probs[i]->compute_constraints(c1,x);
		probs_new[i]->compute_constraints(c2,x);

		if (std::equal(f1.begin(),f1.end(),f2.begin())) {
			std::cout << ": Fitness pass,";
		} else {
			std::cout << ": Fitness FAILED, " << std::endl;
			std::cout << x<< std::endl;
			std::cout << f1 << " " << f2 << std::endl;
			std::cout << *probs[i]<<std::endl;
			std::cout << *probs_new[i]<<std::endl;
			return 1;
		}
		if (std::equal(c1.begin(),c1.end(),c2.begin())) {
			std::cout << " Constraints pass";
		} else {
			std::cout << " Constraints FAILED" << std::endl;
			std::cout << " c1 = " << c1 << std::endl;
			std::cout << " c2 = " << c2 << std::endl;
			return 1;
		}
		}

	}
	std::cout << std::endl;
	return 0;
}
