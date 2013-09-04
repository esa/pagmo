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

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "../src/pagmo.h"

using namespace pagmo;
int main()
{
	unsigned int pop_size = 20;
	unsigned int gen = 123;

	//increment this if you add a multiobjective problem
	unsigned int n_mo = 2;

	// create two containers of pagmo::algorithms
	std::vector<algorithm::base_ptr> algos;
	std::vector<algorithm::base_ptr> algos_new;
	
	// fill it up with algorithm (first the multiobjective ones

	algos.push_back(algorithm::nsga2(gen,0.5,11,0.3,11).clone());
	algos_new.push_back(algorithm::nsga2().clone());
	algos.push_back(algorithm::vega(gen,.9,.021,1,algorithm::vega::mutation::RANDOM,0.3,algorithm::vega::crossover::BINOMIAL).clone());
	algos_new.push_back(algorithm::vega().clone());
	algos.push_back(algorithm::bee_colony(gen,19).clone());
	algos_new.push_back(algorithm::bee_colony().clone());
	algos.push_back(algorithm::cmaes(gen,0.5, 0.5, 0.5, 0.5, 0.7, 1e-5, 1e-5, false).clone());
	algos_new.push_back(algorithm::cmaes().clone());
	algos.push_back(algorithm::cs(gen*10,0.02,0.3,0.3).clone());
	algos_new.push_back(algorithm::cs().clone());
	algos.push_back(algorithm::de(gen,0.9,0.9,3).clone());
	algos_new.push_back(algorithm::de().clone());
	algos.push_back(algorithm::de_1220(1,2,std::vector<int>(1,9),false,1e-5,1e-5).clone());
	algos_new.push_back(algorithm::de_1220().clone());
	algos.push_back(algorithm::ihs(gen,0.2,0.2,0.2,0.2,0.2).clone());
	algos_new.push_back(algorithm::ihs().clone());
	algos.push_back(algorithm::jde(gen,7,2).clone());
	algos_new.push_back(algorithm::jde().clone());
	algos.push_back(algorithm::mbh(algorithm::de(gen),2,0.03).clone());
	algos_new.push_back(algorithm::mbh().clone());
	algos.push_back(algorithm::mde_pbx(gen,0.5,0.5,1e-10,1e-10).clone());
	algos_new.push_back(algorithm::mde_pbx().clone());
	algos.push_back(algorithm::monte_carlo(gen).clone());
	algos_new.push_back(algorithm::monte_carlo().clone());
	algos.push_back(algorithm::ms(algorithm::monte_carlo(gen),5).clone());
	algos_new.push_back(algorithm::ms().clone());
	algos.push_back(algorithm::null().clone());
	algos_new.push_back(algorithm::null().clone());
	algos.push_back(algorithm::pso(gen,0.5,0.5,0.5,0.5,3,3,3).clone());
	algos_new.push_back(algorithm::pso().clone());
	algos.push_back(algorithm::pso_generational(gen,0.5,0.5,0.5,0.5,3,3,3).clone());
	algos_new.push_back(algorithm::pso_generational().clone());
	algos.push_back(algorithm::sa_corana(gen*100,5.0,1e-5,25,10,0.5).clone());
	algos_new.push_back(algorithm::sa_corana().clone());
	algos.push_back(algorithm::sga(gen,.9, .021, 5, algorithm::sga::mutation::RANDOM, 0.3, algorithm::sga::selection::BEST20, algorithm::sga::crossover::BINOMIAL).clone());
	algos_new.push_back(algorithm::sga().clone());

#ifdef PAGMO_ENABLE_GSL
	algos.push_back(algorithm::gsl_bfgs(gen,1e-3,1e-3,0.03,1e-3).clone());
	algos_new.push_back(algorithm::gsl_bfgs().clone());
	algos.push_back(algorithm::gsl_bfgs2(gen,1e-3,1e-3,0.03,1e-3).clone());
	algos_new.push_back(algorithm::gsl_bfgs2().clone());
	algos.push_back(algorithm::gsl_fr(gen,1e-3,1e-3,0.03,1e-3).clone());
	algos_new.push_back(algorithm::gsl_fr().clone());
	algos.push_back(algorithm::gsl_nm(gen,1e-3,0.1).clone());
	algos_new.push_back(algorithm::gsl_nm().clone());
	algos.push_back(algorithm::gsl_nm2(gen,1e-3,0.1).clone());
	algos_new.push_back(algorithm::gsl_nm2().clone()); 
	algos.push_back(algorithm::gsl_nm2rand(gen,1e-3,0.1).clone());
	algos_new.push_back(algorithm::gsl_nm2rand().clone());
	algos.push_back(algorithm::gsl_pr(gen,1e-3,1e-3,0.03,1e-3).clone());
	algos_new.push_back(algorithm::gsl_pr().clone());
#endif
	
#ifdef PAGMO_ENABLE_NLOPT
	algos.push_back(algorithm::nlopt_aug_lag(2, gen, 2.3E-6, 2.3E-6, 150, 2.3E-6, 2.3E-6).clone());
	algos_new.push_back(algorithm::nlopt_aug_lag().clone());
	algos.push_back(algorithm::nlopt_aug_lag_eq(2, gen, 2.3E-6, 2.3E-6, 150, 2.3E-6, 2.3E-6).clone());
	algos_new.push_back(algorithm::nlopt_aug_lag_eq().clone());
	algos.push_back(algorithm::nlopt_bobyqa(gen, 2.3E-6, 2.3E-6).clone());
	algos_new.push_back(algorithm::nlopt_bobyqa().clone());
	algos.push_back(algorithm::nlopt_cobyla(gen, 2.3E-6, 2.3E-6).clone());
	algos_new.push_back(algorithm::nlopt_cobyla().clone());
	algos.push_back(algorithm::nlopt_slsqp(gen, 2.3E-6, 2.3E-6).clone());
	algos_new.push_back(algorithm::nlopt_slsqp().clone());
	algos.push_back(algorithm::nlopt_sbplx(gen, 2.3E-6, 2.3E-6).clone());
	algos_new.push_back(algorithm::nlopt_sbplx().clone());
	algos.push_back(algorithm::nlopt_mma(gen, 2.3E-6, 2.3E-6).clone());
	algos_new.push_back(algorithm::nlopt_mma().clone());
#endif
	
#ifdef PAGMO_ENABLE_SNOPT
	algos.push_back(algorithm::snopt(gen, 2.3E-6, 2.3E-6).clone());
	algos_new.push_back(algorithm::snopt().clone());
#endif

// IPOPT algorithm.
#ifdef PAGMO_ENABLE_IPOPT
	algos.push_back(algorithm::ipopt(gen, 2.3E-6, 2.3E-6, 2.3E-6, false, 2.0, 0.234).clone());
	algos_new.push_back(algorithm::ipopt().clone());
#endif
	// pick a box-constrained, single objective, continuous problem
	problem::ackley prob(10);
	
	// and pick a multiobjective one
	problem::zdt1 prob_mo(10);

	// make a population out of it
	population pop_original(prob,pop_size);
	population pop_original_mo(prob_mo,pop_size);

	//serialize algos and deserialize into algos_new checking they are then identical
	for (size_t i=0; i< algos.size(); ++i) {
		{
		// create and open a character archive for output
		std::ofstream ofs("test.ar");
		// save data to archive
		boost::archive::text_oarchive oa(ofs);
		// write class instance to archive
		oa & algos[i];
		// archive and stream closed when destructors are called
		}
	
		{
		// create and open an archive for input
		std::ifstream ifs("test.ar");
		boost::archive::text_iarchive ia(ifs);
		// read class state from archive
		ia & algos_new[i];
		// archive and stream closed when destructors are called
		}	


		{
		//copy the original population
		population pop1(prob), pop2(prob);
		if (i<n_mo) {
			pop1 = population(pop_original_mo); 
			pop2 = population(pop_original_mo);
		}
		else {
			pop1 = population(pop_original);
			pop2 = population(pop_original);
		}
		algos[i]->evolve(pop1);
		algos_new[i]->evolve(pop2);
		std::cout << std::endl << std::setw(80) << algos[i]->get_name();
		decision_vector x1(pop1.champion().x), x2(pop2.champion().x);
		if (std::equal(x1.begin(),x1.end(),x2.begin())) {
			std::cout << ": pass";
		} else {
			std::cout << ": Champion is different:" << std::endl;
			std::cout << x1 << std::endl;
			std::cout << x2 << std::endl;
			return 1;
		}
		}
	}
	std::cout << std::endl;
	return 0;
}
