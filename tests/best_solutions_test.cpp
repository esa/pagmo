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
#include "test.h"

//--------------------------------------------------------------------------------
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

const double EPS = 1e-6;

///The idea of this unit test is to go through all pagmo::problems, initialize the best solution,
///retrieve the decision vector

using namespace pagmo;

struct best_solutions_test{
	problem::base_ptr problem;
	std::vector<fitness_vector> best_f;
	std::vector<constraint_vector> best_c;
};

best_solutions_test make_struct(problem::base_ptr pb, const std::vector<fitness_vector>& f = std::vector<fitness_vector>(), const std::vector<constraint_vector>& c = std::vector<fitness_vector>()){
	best_solutions_test test;
	test.problem = pb;
	test.best_f = f;
	test.best_c = c;
	return test;
}

std::vector<fitness_vector> get_cec2006_best_f(int id);
std::vector<fitness_vector> get_cec2006_best_c(int id, int c_dimension);

int main()
{
	unsigned int dimension = 24;

	// create a containers of pagmo::problems
	std::vector<best_solutions_test> best_tests;

	// fill it up with problems and best known fitness and cstrs values
	best_tests.push_back(make_struct(problem::ackley(dimension).clone()));
	best_tests.push_back(make_struct(problem::rosenbrock(dimension).clone()));
	best_tests.push_back(make_struct(problem::branin().clone()));
	best_tests.push_back(make_struct(problem::dejong(dimension).clone()));
	best_tests.push_back(make_struct(problem::fon().clone()));
	best_tests.push_back(make_struct(problem::golomb_ruler(10,20).clone()));
	best_tests.push_back(make_struct(problem::griewank(dimension).clone()));
	best_tests.push_back(make_struct(problem::himmelblau().clone()));
	best_tests.push_back(make_struct(problem::string_match("e dai dai dai.....portiamolo a casa!!").clone()));
	best_tests.push_back(make_struct(problem::inventory(7,8,1234).clone()));
	best_tests.push_back(make_struct(problem::kur(dimension).clone()));
	best_tests.push_back(make_struct(problem::lennard_jones(dimension).clone()));
	best_tests.push_back(make_struct(problem::levy5(dimension).clone()));
	best_tests.push_back(make_struct(problem::luksan_vlcek_1(dimension).clone()));
	best_tests.push_back(make_struct(problem::luksan_vlcek_2(dimension).clone()));
	best_tests.push_back(make_struct(problem::luksan_vlcek_3(dimension).clone()));
	best_tests.push_back(make_struct(problem::michalewicz(dimension).clone()));
	best_tests.push_back(make_struct(problem::pol().clone()));
	best_tests.push_back(make_struct(problem::rastrigin(dimension).clone()));
	best_tests.push_back(make_struct(problem::sch().clone()));
	best_tests.push_back(make_struct(problem::schwefel(dimension).clone()));
	best_tests.push_back(make_struct(problem::snopt_toyprob().clone()));

	//----- Test ZDT -----//
	for(int i=1; i<=6;i++) {
	    best_tests.push_back(make_struct(problem::zdt(i, dimension).clone()));
	}

	//----- Test DTLZ -----//
	for(int i=1; i<=7;i++) {
	    best_tests.push_back(make_struct(problem::dtlz(i, dimension).clone()));
	}
	
	//----- Test CEC2006 -----//
	for(int i=1; i<=24; i++){
		best_tests.push_back(make_struct(problem::cec2006(i).clone(), get_cec2006_best_f(i), get_cec2006_best_c(i, problem::cec2006(i).get_c_dimension())));
	}

	//----- Test CEC2009 - UF set-----//
	for(int i=1; i<=10; i++){
		best_tests.push_back(make_struct(problem::cec2009(i, dimension, false).clone()));
	}
	//----- Test CEC2009 - CF set-----//
	for(int i=1; i<=10; i++){
		best_tests.push_back(make_struct(problem::cec2009(i, dimension, true).clone()));
	}

	//----- Test meta-problems -----//
	problem::zdt zdt1_before_transform1(1, dimension);
	//----- shifted -----//
	best_tests.push_back(make_struct(problem::shifted(zdt1_before_transform1).clone()));
	//----- rotated -----//
	best_tests.push_back(make_struct(problem::rotated(zdt1_before_transform1).clone()));

#ifdef PAGMO_ENABLE_KEP_TOOLBOX
	best_tests.push_back(make_struct(problem::cassini_1(2).clone()));
	best_tests.push_back(make_struct(problem::cassini_2().clone()));
	best_tests.push_back(make_struct(problem::gtoc_1().clone()));
	best_tests.push_back(make_struct(problem::messenger().clone()));
	best_tests.push_back(make_struct(problem::rosetta().clone()));
	best_tests.push_back(make_struct(problem::messenger_full().clone()));
	best_tests.push_back(make_struct(problem::tandem(3,10).clone()));
	best_tests.push_back(make_struct(problem::laplace(std::vector<int>(default_sequence,default_sequence + 5)).clone()));
	best_tests.push_back(make_struct(problem::mga_1dsm_alpha(construct_sequence()).clone()));
	best_tests.push_back(make_struct(problem::mga_1dsm_tof(construct_sequence()).clone()));
#endif

	// initialize the best solution and retrieve the decision vector
	for (size_t i=0; i<best_tests.size(); ++i) {
		std::cout << std::endl << std::setw(40) << best_tests[i].problem->get_name() << std::endl;

		if(best_tests[i].problem->get_best_x().empty())
			std::cout << "Best decision vector is not implemented." << std::endl;
		else {
			const std::vector<decision_vector> &best_x = best_tests[i].problem->get_best_x();
			//browse the solutions set
			for(std::vector<decision_vector>::size_type j=0; j<best_x.size(); j++) {
				if(best_tests[i].best_f.empty()) {
					std::cout << "Best fitness vector is not implemented." << std::endl;
				}
				else{
					const fitness_vector &best_f = best_tests[i].best_f.at(j);
					fitness_vector best_f_computed = best_tests[i].problem->get_best_f().at(j);

					if(is_eq_vector(best_f_computed, best_f, EPS)){
						std::cout << " fitness passes, ";
					}
					else{
						std::cout << " fitness failed!"<<std::endl;
						return 1;
					}
				}

				// If the test-local best constraint vector is empty OR 
				// the problem has a constraint vector dimension of zero,
				// then the constraint vector is not implemented and this 
				// check will be skipped.
				if(best_tests[i].best_c.empty() || best_tests[i].problem->get_c_dimension() == 0) {
					std::cout << "Best constraint vector is not implemented." << std::endl;
				}
				else {
					const constraint_vector &best_c = best_tests[i].best_c.at(j);
					constraint_vector best_c_computed = best_tests[i].problem->get_best_c().at(j);

					if(is_eq_vector(best_c_computed, best_c, EPS)){
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


std::vector<fitness_vector> get_cec2006_best_f(int id){
	std::vector<decision_vector> best_f;
	int f_dimension = 1;
	double f = 0.0;
	switch(id)
	{
	case 1:
	{
		f = -15.;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 2:
	{
		f = -0.80361910412559;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 3:
	{
		f = -1.00050010001000;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 4:
	{
		f = -3.066553867178332e+004;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 5:
	{
		f = 5126.4967140071;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 6:
	{
		f = -6961.81387558015;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 7:
	{
		f = 24.30620906818;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 8:
	{
		f = -0.0958250414180359;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 9:
	{
		f = 680.630057374402;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 10:
	{
		f = 7049.24802052867;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 11:
	{
		f = 0.7499;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 12:
	{
		f = -1.;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 13:
	{
		f = 0.053941514041898;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 14:
	{
		f = -47.7648884594915;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 15:
	{
		f = 961.715022289961;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 16:
	{
		f = -1.90515525853479;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 17:
	{
		f = 8853.53967480648;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 18:
	{
		f = -0.866025403784439;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 19:
	{
		f = 32.6555929502463;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 20:
	{
		f = 0.20497940028563599;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 21:
	{
		f = 193.724510070035;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 22:
	{
		f = 236.430975504001;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 23:
	{
		f = -400.055099999999584;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	case 24:
	{
		f = -5.50801327159536;
		fitness_vector fitness(f_dimension,f);
		best_f.push_back(fitness);
		break;
	}
	default:
		pagmo_throw(value_error, "Error: There are only 24 test functions in the CEC2006 test suite!");
		break;
	}

	return best_f;
}

std::vector<constraint_vector> get_cec2006_best_c(int id, int c_dimension){
	std::vector<constraint_vector> best_c;
	switch(id)
	{
	case 1:
	{
		const double c_vector[] = {0., 0., 0., -5., -5., -5., 0., 0., 0.};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 2:
	{
		const double c_vector[] = {-1.2878587085651816e-14, -120.06741615259264};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 3:
	{
		const double c_vector[] = {9.9999999999988987e-05};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 4:
	{
		const double c_vector[] = {0., -92., -11.159499691073137, -8.8405003089268632, -4.9999999999999964, -3.5527136788005009e-15};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 5:
	{
		const double c_vector[] = {9.9999999974897946e-05, 9.9999999974897946e-05, 9.9999999974897946e-05, -0.034890145690411378, -1.0651098543095887};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 6:
	{
		const double c_vector[] = {0., 0.};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 7:
	{
		const double c_vector[] = {5.6843418860808015e-14, -1.1723955140041653e-13, 3.907985046680551e-14, -6.0254023992456496e-12,
								   -7.1054273576010019e-15, -2.8421709430404007e-14, -6.1485036896036398, -50.023961731838071};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 8:
	{
		const double c_vector[] = {-1.737459723297992, -0.16776326380511744};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 9:
	{
		const double c_vector[] = {-1.3766765505351941e-14, -252.56171634346606, -144.87817845461515, -2.4868995751603507e-14};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 10:
	{
		const double c_vector[] = {0, 0, -5.5511151231257827e-16, -1.4551915228366852e-11, -2.9103830456733704e-11, -1.1641532182693481e-10};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 11:
	{
		const double c_vector[] = {9.9999999999988987e-05};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 12:
	{
		const double c_vector[] = {-0.0625};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 13:
	{
		const double c_vector[] = {9.9999999994437871e-05, -0.00010000000000331966, 9.9999999998878764e-05};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 14:
	{
		const double c_vector[] = {0.00010000000000109921, 9.9999999999544897e-05, 0.00010000000000043308};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 15:
	{
		const double c_vector[] = {9.9999999999766942e-05, 9.9999999989108801e-05};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 16:
	{
		const double c_vector[] = {-81.790353283003498, 0, 0, -5.6843418860808015e-14,
								   0, -192.13000000000002, -0.29331594739664624, -1035.8683840526032, -0.094075252933153664,
								   -23.660924747066847, -1.7844298057299284, -449.5725701942701, -337.68448248701077,
								   -239.32051751298917, -131.68954952441663, -133.26545047558335, -1.0081273474015813,
								   -4.4258726525984189, -0.0010113580124478661, -0.074988641987552146, -46.616183234373679,
								   -118.75981676562631, -47.15206111982161, -316.25993888017842, -187.97290376436615,
								   -329.2410962356339, -0.13924565212261086, -518.23575434787733, -1859.6857829062783,
								   -315.19021709372191, -5850.0391214556676, -12032.598878544331, -0.18271094712517674,
								   -0.14028905287482327, -68915.669999999998, 0, -9275494.7446969859, -67900.25530301407};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 17:
	{
		const double c_vector[] = {9.5279025970285147e-05, 9.9999999889632818e-05, -7.5298006409596496e-05, 9.9999999605415724e-05};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 18:
	{
		const double c_vector[] = {0, -0.6402463624140452, 0, 0,
								   -0.64024636445565331, -1.5959455978986625e-16, -7.6327832942979512e-17, -0.64024636141023228,
								   0, -0.67204348836757655, -0.19398191412316787, -0.39453065073841909, -0.47149475433971377};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 19:
	{
		const double c_vector[] = {1.7763568394002505e-15, -0, -0, 8.8817841970012523e-16, -0};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 20:
	{
		const double c_vector[] = {9.9999999999766942e-05, -2.1972568189549347e-17, 1.5539704692828344e-18, 8.1753737821131767e-19,
								   9.9999999999766942e-05, -1.0383082213426958e-17, -2.992491915147666e-17, -3.20853810872953e-18,
								   0, 9.9999999999877964e-05, 9.9999999999840031e-05, 6.8665220979881771e-16,
								   9.9999999999988987e-05, -9.9999999999322853e-05, 0.14375363724895993, 1.7660316449634668e-19,
								   7.5785257203801829e-19, 2.2240466906989736e-18, 2.0929449795679412e-18, 0};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 21:
	{
		const double c_vector[] = {9.9999997473787516e-05, 0.00010000000111176632, 9.9999999999766942e-05, -9.9999999999766942e-05,
							 -9.9999999999766942e-05, -2.8421709430404007e-14};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 22:
	{
		const double c_vector[] = {-2.9260292649269104e-05, 9.1973692178726196e-05, -5.184859037399292e-05, -4.9985945224761963e-05,
								   7.2941184043884277e-05, 9.5665454864501953e-06, 5.6096352636814117e-05, 1.6726553440093994e-05,
								   5.4139643907546997e-05, -7.7500534302998858e-05, 1.6723476818469862e-05, 5.7372450426917965e-05,
								   5.7400810546504033e-05, 7.4797550418281844e-05, 7.4730862352545557e-05, -6.1559895403462406e-07,
								   6.1231858126120642e-06, -8.7799054654169595e-05, 9.5923077310544613e-05, -2.2079848349676467e-07};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 23:
	{
		const double c_vector[] = {-0.00010000000000331966, 9.9999999999988987e-05, -9.9999999999507005e-05, -9.9999999974897946e-05,
								   -2.5000000001256605e-06, 0};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	case 24:
	{
		const double c_vector[] = {-4.8849813083506888e-14, 1.7053025658242404e-13};
		constraint_vector constraint(c_dimension);
		std::copy(c_vector,c_vector + c_dimension,constraint.begin());
		best_c.push_back(constraint);
		break;
	}
	default:
	{
		pagmo_throw(value_error, "Error: There are only 24 test functions in this test suite!");
		break;
	}
	}
	return best_c;
}
